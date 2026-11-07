"""Run two frozen published numerical controls without phenotype prediction."""
import argparse
from decimal import Decimal
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def module(name):
    source = Path(__file__).with_name(name)
    data = source.read_bytes()
    spec = importlib.util.spec_from_file_location(name, source)
    value = importlib.util.module_from_spec(spec)
    exec(compile(data, str(source), "exec"), value.__dict__)
    value.__source_sha256__ = hashlib.sha256(data).hexdigest()
    return value


def printed_comparison(value, expected, places):
    value, expected = Decimal(str(value)), Decimal(expected)
    half_unit = Decimal(1).scaleb(-places) / 2
    lower, upper = expected - half_unit, expected + half_unit
    status = ("published_precision_match" if lower < value < upper else
              "rounding_boundary_ambiguous" if value in (lower, upper) else "discrepancy")
    return {"status": status, "published_similarity_percent": str(expected),
            "printed_decimal_places": places,
            "implied_rounding_interval_percent": {"lower": str(lower), "upper": str(upper),
                                                  "exact_endpoints": "ambiguous; paper rounding convention not established"}}


def run(args):
    policy = json.loads(args.policy.read_text())
    if "validator_sha256" in policy and sha(Path(__file__)) != policy["validator_sha256"]:
        raise ValueError("Validator differs from frozen diagnostic")
    install_provenance = None
    if "tool_install_provenance" in policy:
        install = policy["tool_install_provenance"]
        if sha(Path(install["path"])) != install["sha256"]:
            raise ValueError("Tool installation provenance changed")
        install_provenance = json.loads(Path(install["path"]).read_text())
    previous = None
    if args.reference_report is not None:
        if sha(args.reference_report) != policy["previous_reference_report_sha256"]:
            raise ValueError("Previous reference report differs from frozen version diagnostic")
        previous = json.loads(args.reference_report.read_text())
        if previous["policy_sha256"] != policy["previous_policy_sha256"]:
            raise ValueError("Previous reference policy differs from frozen version diagnostic")
        previous_expectations = {p["id"]: (p["accessions"], p["published_similarity_percent"], p["printed_decimal_places"])
                                 for p in previous["control_results"]}
        if previous_expectations != {p["id"]: (p["accessions"], p["published_similarity_percent"], p["printed_decimal_places"])
                                     for p in policy["control_pairs"]}:
            raise ValueError("Version diagnostic must keep the same reference expectations")
    if sha(args.adapter) != policy["adapter_sha256"]:
        raise ValueError("Adapter differs from frozen reference-check policy")
    if args.workdir.exists() or args.output.exists():
        raise ValueError("Use a new raw work directory and report; preserve prior results")
    roots = [args.policy, args.adapter, args.target_manifest, args.t7_fasta, args.blast_bin]
    for source in roots:
        for destination in (args.workdir, args.output):
            if destination.resolve() == source.resolve() or destination.resolve().is_relative_to(source.resolve()):
                raise ValueError("Outputs must be separate from inputs")
    normalizer, blast_helper = module("normalize-genophi-phages.py"), module("audit-genophi-phage-overlap.py")
    for name, imported in (("normalizer_sha256", normalizer), ("completion_helper_sha256", blast_helper)):
        if name in policy and imported.__source_sha256__ != policy[name]:
            raise ValueError("Validation helper differs from frozen diagnostic")
    args.workdir.mkdir(parents=True)
    shutil.copyfile(args.policy, args.workdir / "frozen-policy.json")
    adapter = args.workdir / "frozen-adapter.py"
    shutil.copyfile(args.adapter, adapter)
    commands, tool_info = [], {}
    for name in ("blastn", "makeblastdb"):
        executable = args.blast_bin / name
        version = subprocess.check_output([str(executable), "-version"], text=True).splitlines()[0]
        if version.split(": ")[-1] != policy["local_blast_version"]:
            raise ValueError(f"Unexpected BLAST version: {version}")
        binary_sha = sha(executable)
        if "executable_sha256" in policy and binary_sha != policy["executable_sha256"][name]:
            raise ValueError("BLAST executable differs from frozen diagnostic")
        tool_info[name] = {"version": version, "sha256": binary_sha, "path": str(executable)}
    target_spec = json.loads(args.target_manifest.read_text())["genome_sources"]["T7"]
    sequences, sources, failures = {}, {}, {}
    for accession in sorted({a for pair in policy["control_pairs"] for a in pair["accessions"]}):
        path = args.workdir / "references" / f"{accession}.fna"
        path.parent.mkdir(exist_ok=True)
        attempts = []
        try:
            if previous is not None:
                original = previous["sources"][accession]
                frozen_source = Path(original["path"])
                if sha(frozen_source) != original["file_sha256"]:
                    raise ValueError("Previously frozen reference FASTA changed")
                shutil.copyfile(frozen_source, path)
                source_url = original["url"]
                attempts.append({"source": str(frozen_source), "operation": "copy previously frozen public reference"})
            elif accession == "NC_001604.1":
                shutil.copyfile(args.t7_fasta, path)
                source_url = target_spec["url"]
                attempts.append({"source": str(args.t7_fasta), "operation": "copy existing public reference"})
            else:
                source_url = (f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={accession}&rettype=fasta&retmode=text"
                              if accession.startswith("NC_") else f"https://www.ebi.ac.uk/ena/browser/api/fasta/{accession}?download=false")
                temporary = path.with_suffix(".part")
                command = ["curl", "--fail", "--silent", "--show-error", "--location", "--max-time", "60",
                           "--retry", "1", "--retry-delay", "1", "--max-filesize", "2000000",
                           source_url, "--output", str(temporary)]
                commands.append(command)
                result = subprocess.run(command, capture_output=True, text=True, timeout=150)
                attempts.append({"url": source_url, "returncode": result.returncode, "stderr": result.stderr.strip()})
                if result.returncode:
                    temporary.unlink(missing_ok=True)
                    raise ValueError("Public FASTA retrieval failed; see recorded attempt")
                temporary.replace(path)
            data = path.read_bytes()
            header = data.decode("ascii").splitlines()[0]
            if not header.startswith(">") or accession not in re.split(r"[|\s]+", header[1:]):
                raise ValueError("FASTA header does not establish the exact accession version")
            sequence, changes = normalizer.parse_complete_genome(data)
            sequence_sha = hashlib.sha256(sequence.encode()).hexdigest()
            if previous is not None and (sequence_sha != original["sequence_sha256"] or len(sequence) != original["length"]):
                raise ValueError("Version diagnostic reference differs from previous sequence")
            if accession == "NC_001604.1" and (sequence_sha != target_spec["sequence_sha256"] or len(sequence) != target_spec["length"]):
                raise ValueError("T7 sequence differs from frozen WIsH reference")
            sequences[accession] = sequence
            sources[accession] = {"accession": accession, "url": source_url, "header": header,
                                  "path": str(path), "file_sha256": sha(path), "sequence_sha256": sequence_sha,
                                  "length": len(sequence), "attempts": attempts, **changes}
            if accession == "NC_001604.1":
                sources[accession]["wish_manifest_sequence_verified"] = True
        except (ValueError, OSError, subprocess.SubprocessError) as exc:
            path.with_suffix(".part").unlink(missing_ok=True)
            failures[accession] = {"error": str(exc), "attempts": attempts}
    input_manifest = {"sources": sources, "failures": failures, "policy_sha256": sha(args.policy)}
    (args.workdir / "input-manifest.json").write_text(json.dumps(input_manifest, indent=2) + "\n")
    pairs = []
    for expected in policy["control_pairs"]:
        record = {key: expected[key] for key in ("id", "accessions", "published_similarity_percent", "printed_decimal_places", "locator")}
        if any(a not in sequences for a in expected["accessions"]):
            pairs.append({**record, "status": "blocked", "reason": "One or more exact-version source inputs failed retrieval or authentication"})
            continue
        directory = args.workdir / expected["id"]
        directory.mkdir()
        fasta, lengths, raw, numeric, output = [directory / name for name in
            ("genomes.fna", "lengths.json", "blastn.tsv", "blastn-data.tsv", "similarities.json")]
        accessions = expected["accessions"]
        fasta.write_text("".join(f">r{i}\n{sequences[a]}\n" for i, a in enumerate(accessions)))
        lengths.write_text(json.dumps({f"r{i}": len(sequences[a]) for i, a in enumerate(accessions)}, indent=2) + "\n")
        database = directory / "db"
        pair_commands = [
            [str(args.blast_bin / "makeblastdb"), "-in", str(fasta), "-dbtype", "nucl", "-parse_seqids", "-out", str(database)],
            [str(args.blast_bin / "blastn"), "-query", str(fasta), "-db", str(database), *policy["blastn_parameters"],
             "-outfmt", policy["blastn_format"], "-out", str(raw)],
        ]
        try:
            for index, command in enumerate(pair_commands):
                commands.append(command)
                result = subprocess.run(command, capture_output=True, text=True, timeout=300)
                (directory / f"step-{index}.log").write_text(result.stdout + result.stderr)
                if result.returncode:
                    raise ValueError(f"Reference step {index} failed with exit status {result.returncode}")
            rows, completion = blast_helper.completed_blast(raw.read_text(), {"r0", "r1"}, 14)
            if {row[0] for row in rows if row[0] == row[1]} != {"r0", "r1"}:
                raise ValueError("Reference self-hit controls incomplete")
            numeric.write_text("".join("\t".join(row) + "\n" for row in rows))
            command = [sys.executable, str(adapter), "--blast-tsv", str(numeric), "--genome-lengths", str(lengths), "--output", str(output)]
            commands.append(command)
            result = subprocess.run(command, capture_output=True, text=True, timeout=60)
            (directory / "step-2.log").write_text(result.stdout + result.stderr)
            if result.returncode:
                raise ValueError("Standalone adaptation failed; inspect step-2.log")
            calculated = json.loads(output.read_text())
            if calculated["upstream_equivalence_verified"] is not False or len(calculated["pairs"]) != 1:
                raise ValueError("Unexpected standalone adaptation output")
            pair = calculated["pairs"][0]
            comparison = printed_comparison(pair["similarity_percent_unrounded"], expected["published_similarity_percent"], expected["printed_decimal_places"])
            pairs.append({**record, **comparison, "calculated": pair, "query_completion": completion,
                          "raw_directory": str(directory), "artifacts": {name: sha(directory / name) for name in
                              ("genomes.fna", "lengths.json", "blastn.tsv", "blastn-data.tsv", "similarities.json")}})
        except (ValueError, OSError, subprocess.SubprocessError) as exc:
            pairs.append({**record, "status": "blocked", "reason": str(exc), "raw_directory": str(directory)})
    (args.workdir / "commands.json").write_text(json.dumps(commands, indent=2) + "\n")
    completed = sum(p["status"] != "blocked" for p in pairs)
    return {"schema_version": 1, "purpose": policy["purpose"], "primary_source": policy["primary_source"],
            "policy_sha256": sha(args.policy), "validator_sha256": sha(Path(__file__)),
            "adapter_sha256": sha(adapter), "completion_helper_sha256": blast_helper.__source_sha256__,
            "normalizer_sha256": normalizer.__source_sha256__, "tool_install_provenance": install_provenance,
            "wish_manifest_sha256": sha(args.target_manifest), "tools": tool_info,
            "input_manifest_sha256": sha(args.workdir / "input-manifest.json"), "sources": sources,
            "source_failures": failures, "control_results": pairs, "completed_controls": completed,
            "all_controls_match_published_precision": completed == len(pairs) and all(p["status"] == "published_precision_match" for p in pairs),
            "commands": commands, "commands_sha256": sha(args.workdir / "commands.json"),
            "paper_blast_version": policy["paper_blast_version"], "local_blast_version": policy["local_blast_version"],
            "previous_reference_report_sha256": sha(args.reference_report) if previous is not None else None,
            "previous_2_17_results_retained": previous is not None,
            "previous_control_results": [{"id": p["id"], "status": p["status"], "calculated": p.get("calculated")}
                                          for p in previous["control_results"]] if previous is not None else [],
            "parameters_retuned": False, "upstream_equivalence_verified": False,
            "phenotype_labels_loaded": False, "infection_or_assembly_predictions_produced": False,
            "interpretation": policy["interpretation"]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ("policy", "blast-bin", "adapter", "target-manifest", "t7-fasta", "workdir", "output"):
        parser.add_argument("--" + flag, type=Path, required=True)
    parser.add_argument("--reference-report", type=Path,
                        help="Copy authenticated source FASTA from a prior frozen control run; no new sequence retrieval")
    args = parser.parse_args()
    report = run(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    raise SystemExit(0 if report["completed_controls"] == len(report["control_results"]) else 1)
