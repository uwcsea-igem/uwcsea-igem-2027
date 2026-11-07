"""Authenticate and normalize published complete-phage FASTA inputs for audit.

No genome construction: all retained bases must be present in one source record.
The one reviewed empty-record exception is bound to exact source bytes.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import re
from itertools import product

REVIEWED_EMPTY_RECORDS = {
    "3cbcc8a322d3d80a17df4c24a90fd948936005b7bb9ff7effa640fde2ecceab5": [0],
}
REFERENCE_PROTEINS = {"T7": "NP_042005.1", "T3": "NP_523342.1"}


def extract_gp17(data: bytes, accession: str, protein_id: str, genome_sha256: str) -> tuple[str, dict]:
    text = data.decode("ascii")
    version = re.search(r"^VERSION\s+(\S+)", text, re.MULTILINE)
    if version is None or version[1] != accession:
        raise ValueError("GenBank accession differs from frozen reference")
    origin = text.split("\nORIGIN", 1)[1].split("//", 1)[0]
    genome = re.sub(r"[\s\d]", "", origin).upper()
    if hashlib.sha256(genome.encode()).hexdigest() != genome_sha256:
        raise ValueError("GenBank DNA differs from frozen reference sequence")
    features = re.split(r"\n     (?=\S)", text.split("FEATURES", 1)[1].split("\nORIGIN", 1)[0])
    matches = [f for f in features if f.startswith("CDS ") and f'/protein_id="{protein_id}"' in f]
    if len(matches) != 1 or '/product="tail fiber protein"' not in matches[0]:
        raise ValueError("Expected unique annotated tail-fiber CDS")
    feature = matches[0]
    location = re.match(r"CDS\s+(\d+)\.\.(\d+)\n", feature)
    translation = re.search(r'/translation="([^"]+)"', feature)
    if location is None or translation is None:
        raise ValueError("Unsupported reference CDS representation")
    protein = "".join(translation[1].split())
    if not set(protein) <= set("ACDEFGHIKLMNPQRSTVWY"):
        raise ValueError("Ambiguous reference protein")
    start, end = int(location[1]), int(location[2])
    if end - start + 1 != 3 * (len(protein) + 1) or not 1 <= start <= end <= len(genome):
        raise ValueError("Reference CDS/protein lengths disagree")
    # These two reference CDSs have standard ATG initiation, codon_start=1 and
    # bacterial code 11 (standard internal codon meanings). Reject other forms.
    if "/codon_start=" in feature and "/codon_start=1" not in feature:
        raise ValueError("Unsupported codon start")
    if "/transl_table=" in feature and "/transl_table=11" not in feature:
        raise ValueError("Unsupported genetic code")
    codons = ("".join(parts) for parts in product("TCAG", repeat=3))
    residues = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    table = dict(zip(codons, residues))
    cds = genome[start - 1:end]
    if not cds.startswith("ATG") or not set(cds) <= set("ACGT"):
        raise ValueError("Unsupported reference CDS bases/initiation")
    translated = "".join(table[cds[i:i+3]] for i in range(0, len(cds), 3))
    if translated != protein + "*":
        raise ValueError("Annotated reference protein disagrees with CDS translation")
    return protein, {"accession": accession, "protein_id": protein_id, "cds_location": [start, end],
                     "amino_acids": len(protein), "genbank_sha256": hashlib.sha256(data).hexdigest(),
                     "cds_translation_verified": True,
                     "genome_sequence_sha256": genome_sha256,
                     "protein_sequence_sha256": hashlib.sha256(protein.encode()).hexdigest()}


def parse_complete_genome(data: bytes, allowed_empty_records: list[int] | None = None) -> tuple[str, dict]:
    records, current = [], None
    for line in data.decode("ascii").splitlines():
        if not line.strip():
            continue
        if line.startswith(">"):
            if not line[1:].strip():
                raise ValueError("Empty FASTA header")
            records.append([])
            current = records[-1]
        elif current is None:
            raise ValueError("Sequence before FASTA header")
        else:
            sequence = "".join(line.split()).upper()
            if not set(sequence) <= set("ACGTRYSWKMBDHVN"):
                raise ValueError("Non-IUPAC DNA input")
            current.append(sequence)
    sequences = ["".join(parts) for parts in records]
    empty = [i for i, seq in enumerate(sequences) if not seq]
    if empty != (allowed_empty_records or []):
        raise ValueError("Empty-record pattern differs from reviewed source exception")
    retained = [seq for seq in sequences if seq]
    if len(retained) != 1:
        raise ValueError("Require one complete source sequence; never concatenate real contigs")
    return retained[0], {"source_records": len(records), "removed_empty_record_indices": empty,
                         "retained_records": 1, "retained_bases": len(retained[0])}


def run(metadata: Path, genomes: Path, target_genomes: Path, target_manifest: Path, output: Path,
        genbank_references: Path | None = None) -> dict:
    if output.exists():
        raise ValueError("Normalization output must be a new directory")
    spec = importlib.util.spec_from_file_location("prepare", Path(__file__).with_name("prepare-genophi-genomes.py"))
    prepare = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(prepare)
    inventory, provenance = prepare.load_inventory(metadata)
    sources = sorted([item for item in inventory if item["kind"] == "phage"], key=lambda x: x["source_id"])
    normalized, manifest = [], []
    for i, item in enumerate(sources):
        path = genomes / item["path"]
        data = path.read_bytes()
        blob = hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()
        if len(data) != item["expected_bytes"] or blob != item["expected_git_blob_sha1"]:
            raise ValueError(f"Pinned source bytes differ: {item['source_id']}")
        digest = hashlib.sha256(data).hexdigest()
        sequence, changes = parse_complete_genome(data, REVIEWED_EMPTY_RECORDS.get(digest))
        identifier = f"s{i:03d}"
        normalized.append((identifier, sequence))
        manifest.append({"id": identifier, "source_id": item["source_id"], "role": "external_collection",
                         "source_file_sha256": digest, "git_blob_sha1": blob, "url": item["url"],
                         "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(), **changes})
    target_spec = json.loads(target_manifest.read_text())
    for i, name in enumerate(target_spec["phages"]):
        path = target_genomes / "phages" / f"{name}.fna"
        data = path.read_bytes()
        sequence, changes = parse_complete_genome(data)
        digest = hashlib.sha256(sequence.encode()).hexdigest()
        if digest != target_spec["genome_sources"][name]["sequence_sha256"]:
            raise ValueError(f"Target sequence differs from frozen reference: {name}")
        identifier = f"t{i:03d}"
        normalized.append((identifier, sequence))
        manifest.append({"id": identifier, "source_id": name, "role": "target_reference",
                         "accession": target_spec["phages"][name],
                         "source_file_sha256": hashlib.sha256(data).hexdigest(),
                         "sequence_sha256": digest, **changes})
    if len(sources) != 97 or len(normalized) != 99:
        raise ValueError("Expected 97 source genomes and two targets")
    output.mkdir(parents=True)
    fasta = output / "genomes.fna"
    fasta.write_text("".join(f">{identifier}\n{sequence}\n" for identifier, sequence in normalized))
    lengths = output / "lengths.json"
    lengths.write_text(json.dumps({identifier: len(seq) for identifier, seq in normalized}, indent=2) + "\n")
    result = {"schema_version": 1, "purpose": "Public phage sequence-overlap audit inputs; no infection or assembly predictions",
              "normalizer_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "source_provenance": provenance, "target_manifest_sha256": hashlib.sha256(target_manifest.read_bytes()).hexdigest(),
              "genomes_fasta_sha256": hashlib.sha256(fasta.read_bytes()).hexdigest(),
              "lengths_sha256": hashlib.sha256(lengths.read_bytes()).hexdigest(),
              "normalization": "ASCII uppercase, neutral headers, wrapping removed; one exact-hash-reviewed empty record removed; no retained bases changed or contigs joined",
              "genomes": manifest, "external_genomes": len(sources), "target_references": 2,
              "phenotype_labels_loaded": False, "new_sequences_designed": False}
    if genbank_references is not None:
        proteins, protein_manifest = [], {}
        for i, name in enumerate(target_spec["phages"]):
            data = (genbank_references / f"{name}.gb").read_bytes()
            protein, record = extract_gp17(data, target_spec["phages"][name], REFERENCE_PROTEINS[name],
                                          target_spec["genome_sources"][name]["sequence_sha256"])
            proteins.append((f"t{i:03d}_gp17", protein))
            protein_manifest[name] = record
        protein_path = output / "reference-gp17.faa"
        protein_path.write_text("".join(f">{name}\n{seq}\n" for name, seq in proteins))
        result["reference_proteins"] = protein_manifest
        result["reference_gp17_fasta_sha256"] = hashlib.sha256(protein_path.read_bytes()).hexdigest()
    (output / "manifest.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ("metadata-cache", "genomes", "target-genomes", "target-manifest", "output-dir"):
        parser.add_argument("--" + flag, type=Path, required=True)
    parser.add_argument("--genbank-references", type=Path)
    args = parser.parse_args()
    for source in (args.metadata_cache, args.genomes, args.target_genomes):
        if args.output_dir.resolve().is_relative_to(source.resolve()):
            parser.error("Output must be separate from raw inputs")
    run(args.metadata_cache, args.genomes, args.target_genomes, args.target_manifest, args.output_dir, args.genbank_references)
