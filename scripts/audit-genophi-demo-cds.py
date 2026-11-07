"""Translate published demo coordinates against frozen deposited DNA.

No gene caller, model, phenotype table, coordinate search or sequence editing
is involved. The earlier annotation comparison remains unchanged.
"""
import argparse
from collections import Counter
import hashlib
import importlib.util
from itertools import product
import json
from pathlib import Path
import re


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def module(name: str):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).with_name(name))
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def checked(item: dict) -> bytes:
    data = Path(item["path"]).read_bytes()
    if sha(data) != item["sha256"]:
        raise ValueError("Frozen source changed: " + item["path"])
    return data


def parse_genes(data: bytes) -> list[dict]:
    entries, header, chunks = [], None, []
    for line in data.decode("ascii").splitlines() + [">"]:
        if not line.strip():
            continue
        if line.startswith(">"):
            if header is not None:
                match = re.fullmatch(r"(\S+) # (\d+) # (\d+) # (-?1) # (.+)", header)
                if match is None:
                    raise ValueError("Unsupported coordinate header")
                fields = dict(part.split("=", 1) for part in match[5].split(";") if "=" in part)
                sequence = "".join(chunks).upper()
                if not sequence or not set(sequence.removesuffix("*")) <= set("ACDEFGHIKLMNPQRSTVWY"):
                    raise ValueError("Empty/ambiguous protein or internal stop")
                entries.append({"id": match[1], "start": int(match[2]), "end": int(match[3]),
                                "strand": int(match[4]), "partial": fields["partial"],
                                "start_type": fields["start_type"], "code": int(fields["genetic_code"]),
                                "protein": sequence})
            header, chunks = line[1:], []
        elif header is None:
            raise ValueError("Protein before FASTA header")
        else:
            chunks.append("".join(line.split()))
    if not entries or len({g["id"] for g in entries}) != len(entries):
        raise ValueError("Empty proteome or duplicate protein identifiers")
    return entries


def translate_check(genome: str, gene: dict) -> dict:
    if gene["code"] != 11 or gene["partial"] not in {"00", "01", "10", "11"}:
        raise ValueError("Unsupported code or partial flags")
    if not 1 <= gene["start"] <= gene["end"] <= len(genome) or gene["strand"] not in {-1, 1}:
        raise ValueError("Invalid CDS coordinates/strand")
    cds = genome[gene["start"] - 1:gene["end"]]
    if not set(cds) <= set("ACGT") or len(cds) % 3:
        raise ValueError("Ambiguous CDS or incomplete codon")
    if gene["strand"] == -1:
        cds = cds.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    table = dict(zip(("".join(c) for c in product("TCAG", repeat=3)),
                     "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"))
    translated = "".join(table[cds[i:i + 3]] for i in range(0, len(cds), 3))
    partial_start, partial_stop = gene["partial"] if gene["strand"] == 1 else gene["partial"][::-1]
    start_type = gene["start_type"]
    if start_type == "Edge":
        start_consistent = partial_start == "1"
    elif start_type in {"ATG", "GTG", "TTG"}:
        start_consistent = partial_start == "0" and cds[:3] == start_type
        # Code11 initiation has the same M meaning for these three start codons.
        translated = "M" + translated[1:]
    else:
        raise ValueError("Unsupported start type")
    stop_consistent = translated.endswith("*") == (partial_stop == "0")
    no_internal_stop = "*" not in translated.removesuffix("*")
    matches = translated == gene["protein"]
    return {"id": gene["id"], "translation_matches": matches,
            "declared_start_consistent": start_consistent, "declared_stop_consistent": stop_consistent,
            "no_internal_stop": no_internal_stop,
            "translated_protein_sha256": sha(translated.encode()),
            "deposited_protein_sha256": sha(gene["protein"].encode()),
            "passed": matches and start_consistent and stop_consistent and no_internal_stop}


def coordinates(gene: dict) -> tuple[int, int, int]:
    return gene["start"], gene["end"], gene["strand"]


def stop_key(gene: dict) -> tuple[int, int]:
    return gene["strand"], gene["end"] if gene["strand"] == 1 else gene["start"]


def coordinate_comparison(deposited: list[dict], predicted: list[dict]) -> dict:
    left, right = Counter(map(coordinates, deposited)), Counter(map(coordinates, predicted))
    changed = []
    for direction, source, destination in (("deposited", deposited, predicted), ("predicted", predicted, deposited)):
        dest_coords = {coordinates(g) for g in destination}
        dest_stops = {stop_key(g) for g in destination}
        for gene in source:
            if coordinates(gene) not in dest_coords:
                changed.append({"side": direction, "id": gene["id"], "coordinates": list(coordinates(gene)),
                                "same_strand_stop_in_other": stop_key(gene) in dest_stops})
    return {"exact_shared_coordinate_copies": sum((left & right).values()),
            "unshared_coordinates": changed}


def run(policy_path: Path) -> dict:
    policy = json.loads(policy_path.read_bytes())
    inputs = {key: json.loads(checked(value)) for key, value in policy["inputs"].items()}
    helper = module("check-genophi-annotation.py")
    annotation_policy = inputs["annotation_policy"]
    manifest_spec = annotation_policy["phage_manifest"]
    manifest = json.loads(checked(manifest_spec))
    sequence_map = dict(helper.records(checked({"path": manifest_spec["normalized_fasta"],
                                               "sha256": manifest["genomes_fasta_sha256"]})))
    source_map = {g["source_id"]: g for g in manifest["genomes"] if g["role"] == "external_collection"}
    upstream = annotation_policy["upstream_repository"]
    tree = json.loads(Path(upstream["tree_cache"]).read_bytes())
    entries = module("audit-genophi-sources.py").validate_tree(tree, upstream["root_tree_sha1"])
    lengths = {f["source_id"]: f for f in inputs["source_lengths"]["files"]}
    prior = {f["source_id"]: f for f in inputs["annotation_report"]["files"]}
    reports = []
    for example in inputs["annotation_report"]["demo_comparisons"]:
        name, path = example["source_id"], example["path"]
        data = (Path(annotation_policy["demo_comparison"]["directory"]) / Path(path).name).read_bytes()
        blob = hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()
        if blob != entries[path]["sha"] or sha(data) != example["file_sha256"]:
            raise ValueError("Deposited example changed")
        source = source_map[name]
        genome = sequence_map[source["id"]]
        if sha(genome.encode()) != source["sequence_sha256"] or lengths[name]["sequence_sha256"] != source["sequence_sha256"]:
            raise ValueError("Deposited DNA/source-length correspondence differs")
        deposited = parse_genes(data)
        generated_data = checked({"path": str(Path(policy["annotation_workdir"]) / prior[name]["id"] / "proteins.faa"),
                                  "sha256": prior[name]["faa_sha256"]})
        predicted = parse_genes(generated_data)
        results = []
        for gene in deposited:
            try:
                results.append(translate_check(genome, gene))
            except ValueError as error:
                results.append({"id": gene["id"], "passed": False, "error": str(error)})
        metadata_length = lengths[name]["metadata_genome_size"]
        reports.append({"source_id": name, "deposited_faa_git_blob_sha1": blob,
                        "deposited_faa_sha256": sha(data), "source_sequence_sha256": source["sequence_sha256"],
                        "regenerated_faa_sha256": sha(generated_data), "source_bases": len(genome),
                        "metadata_genome_size": metadata_length,
                        "deposited_proteins": len(deposited), "translation_checks_passed": sum(r["passed"] for r in results),
                        "failed_checks": [r for r in results if not r["passed"]],
                        "deposited_coordinates_beyond_metadata_length": sum(g["end"] > metadata_length for g in deposited) if metadata_length else None,
                        "predicted_coordinates_beyond_metadata_length": sum(g["end"] > metadata_length for g in predicted) if metadata_length else None,
                        **coordinate_comparison(deposited, predicted)})
    if len(reports) != 25 or sum(r["deposited_proteins"] for r in reports) != 2980:
        raise ValueError("Require all25 examples and2980 proteins")
    return {"schema_version": 1, "purpose": policy["purpose"], "policy_sha256": sha(policy_path.read_bytes()),
            "helper_sha256": sha(Path(__file__).read_bytes()), "inputs": policy["inputs"],
            "proteomes": len(reports), "deposited_proteins": 2980,
            "translation_checks_passed": sum(r["translation_checks_passed"] for r in reports),
            "coordinate_changes": dict(Counter((g["side"] + ("_same_stop_changed_start" if g["same_strand_stop_in_other"] else "_no_same_stop"))
                                               for r in reports for g in r["unshared_coordinates"])),
            "exact_shared_coordinate_copies": sum(r["exact_shared_coordinate_copies"] for r in reports),
            "files": reports, "gene_callers_rerun": False, "sequences_changed": False,
            "phenotype_labels_loaded": False, "interaction_scores_generated": False,
            "interpretation": policy["interpretation"]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--policy", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("Use a new report path; preserve existing evidence")
    result = run(args.policy)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    raise SystemExit(0 if result["translation_checks_passed"] == result["deposited_proteins"] else 1)
