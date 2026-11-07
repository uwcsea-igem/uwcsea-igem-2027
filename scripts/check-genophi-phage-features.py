"""Frozen, label-free GenoPHI phage representation pilot (no interaction model).

Standard-library adaptation of pinned GenoPHI family/assignment bookkeeping.
MMseqs18 is an explicit software adaptation; raw alignments remain replayable.
"""
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import resource
import signal
import struct
import subprocess
import sys
import time


def digest(data):
    return hashlib.sha256(data).hexdigest()


def checked(path, expected):
    data = Path(path).read_bytes()
    if digest(data) != expected:
        raise ValueError(f"Checksum differs: {path}")
    return data


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def proteins(data):
    """Keep input order and copies; remove exactly one terminal stop."""
    out, name, chunks = [], None, []
    for line in data.decode("ascii").splitlines():
        if not line.strip():
            continue
        if line.startswith(">"):
            if name is not None:
                out.append((name, "".join(chunks).removesuffix("*")))
            if not line[1:].strip():
                raise ValueError("Empty FASTA header")
            name, chunks = line[1:].split()[0], []
        else:
            if name is None or any(c.isspace() for c in line):
                raise ValueError("Malformed protein FASTA")
            chunks.append(line)
    if name is not None:
        out.append((name, "".join(chunks).removesuffix("*")))
    if not out or len({n for n, _ in out}) != len(out):
        raise ValueError("Empty FASTA or duplicate protein ID")
    if any(not seq or not set(seq) <= set("ACDEFGHIKLMNPQRSTVWYBXZJUO") for _, seq in out):
        raise ValueError("Empty protein, internal stop or invalid amino acid")
    return out


def f32(text):
    try:
        value = struct.unpack("f", struct.pack("f", float(text)))[0]
    except (ValueError, OverflowError, struct.error) as exc:
        raise ValueError("Invalid float32 alignment value") from exc
    if not math.isfinite(value):
        raise ValueError("Nonfinite alignment value")
    return value


def read_families(path, train_ids):
    mapping = {}
    for line in Path(path).read_text().splitlines():
        fields = line.split("\t")
        if len(fields) != 2:
            raise ValueError("Cluster TSV must have two columns")
        family, member = fields
        if family not in train_ids or member not in train_ids or member in mapping:
            raise ValueError("Invalid/duplicate/nontraining cluster member")
        mapping[member] = family
    if set(mapping) != set(train_ids) or not set(mapping.values()) <= set(mapping):
        raise ValueError("Incomplete cluster member coverage")
    return mapping


def best_hits(path, queries, train_ids, families):
    """Pinned float32 identity/E-value ordering; stable first raw-row ties."""
    best, tied, rows, seen = {}, {}, 0, Counter()
    for lineno, line in enumerate(Path(path).read_text().splitlines(), 1):
        fields = line.split("\t")
        if len(fields) != 11:
            raise ValueError(f"Expected 11 alignment fields at row {lineno}")
        q, target = fields[:2]
        if q not in queries or target not in train_ids:
            raise ValueError("Alignment contains nonquery or nontraining target")
        score, identity, evalue = map(f32, fields[2:5])
        if not 0 <= identity <= 1 or evalue < 0:
            raise ValueError("Invalid identity/E-value")
        try:
            qs, qe, qlen, ts, te, tlen = map(int, fields[5:])
        except ValueError as exc:
            raise ValueError("Invalid alignment coordinates") from exc
        if not (qlen == len(queries[q]) and tlen == len(train_ids[target])
                and 0 <= qs <= qe < qlen and 0 <= ts <= te < tlen):
            raise ValueError("Alignment coordinates/length differ from proteins")
        key = (-identity, evalue)
        row = {"target": target, "family": families[target], "identity_float32": identity,
               "evalue_float32": evalue, "raw_row": lineno}
        if q not in best or key < best[q][0]:
            best[q] = (key, row)
            tied[q] = {families[target]}
        elif key == best[q][0]:
            tied[q].add(families[target])
        rows += 1
        seen[q] += 1
    ambiguities = [{"query": q, "families": sorted(fs), "chosen": best[q][1]}
                   for q, fs in sorted(tied.items()) if len(fs) > 1]
    return {q: row for q, (_, row) in best.items()}, {
        "raw_rows": rows, "queries_with_hits": len(best),
        "cross_family_best_hit_ties": ambiguities,
        "raw_schema": ["Query", "Target", "Score", "SeqIdentity", "E-value",
                       "qStartPos", "qEndPos", "qLen", "tStartPos", "tEndPos", "tLen"]}


def consolidate(train_names, train_protein_genome, hits):
    occurrence = defaultdict(set)
    for query, hit in hits.items():
        occurrence[hit["family"]].add(train_protein_genome[query])
    retained = {family: tuple(int(g in genomes) for g in train_names)
                for family, genomes in sorted(occurrence.items()) if len(genomes) > 1}
    pattern_features, family_features = {}, {}
    for family, pattern in retained.items():
        if pattern not in pattern_features:
            pattern_features[pattern] = f"pc_{len(pattern_features)}"
        family_features[family] = pattern_features[pattern]
    return family_features, retained, {
        "families_after_self_search": len(occurrence),
        "families_retained_more_than_one_training_genome": len(retained),
        "features_after_exact_occurrence_consolidation": len(pattern_features)}


def project(genomes, protein_genome, hits, family_features):
    features = sorted(set(family_features.values()))
    found = {g: set() for g in genomes}
    for protein, hit in hits.items():
        feature = family_features.get(hit["family"])
        if feature is not None:
            found[protein_genome[protein]].add(feature)
    return features, {g: [int(f in found[g]) for f in features] for g in genomes}


def prepare(policy):
    sources = {item["path"]: json.loads(checked(item["path"], item["sha256"]))
               for item in policy["inputs"]}
    mapping = sources["data/processed/genophi-training-source-map.json"]
    names = policy["training_phages"]
    if (len(names) != 93 or names != sorted(set(names))
            or sorted(set(mapping["phages"]) - {"427_P2"}) != names
            or "T7_Portugal" in mapping["phages"]
            or policy["target_phages"] != ["T7", "T3"]):
        raise ValueError("Frozen train/target membership differs")
    checked(policy["mmseqs_executable"], policy["mmseqs_sha256"])
    annotation = sources["data/processed/genophi-annotation-check.json"]
    aa, verified = {}, []
    for item in annotation["files"]:
        path = Path(policy["annotation_workdir"]) / item["id"] / "proteins.faa"
        data = checked(path, item["faa_sha256"])
        verified.append({"path": str(path), "sha256": item["faa_sha256"]})
        if item["kind"] == "phage":
            if item["source_id"] in aa:
                raise ValueError("Duplicate annotated phage")
            records = proteins(data)
            if len(records) != item["proteins"]:
                raise ValueError("Protein count differs from annotation audit")
            aa[item["source_id"]] = records
    required = names + policy["target_phages"]
    if any(g not in aa for g in required):
        raise ValueError("Missing annotated proteome")
    return aa, verified


def normalized_inputs(policy, aa, workdir):
    records, train, target = [], {}, {}
    train_genome, target_genome, genome_records = {}, {}, []
    selections = [(g, "train", aa[g]) for g in policy["training_phages"]]
    selections += [(g, "target", aa[g]) for g in policy["target_phages"]]
    selections += [("mapping_control", "control", aa[policy["training_phages"][0]])]
    for gi, (genome, role, source) in enumerate(selections):
        normalized = []
        for pi, (original, seq) in enumerate(source):
            neutral = f"g{gi:03d}p{pi:06d}"
            mapping = {"protein": neutral, "genome": genome, "role": role,
                       "source_protein": original, "sequence_sha256": digest(seq.encode()),
                       "length": len(seq)}
            records.append(mapping)
            normalized.append(f">{neutral}\n{seq}\n")
            if role == "train":
                train[neutral], train_genome[neutral] = seq, genome
            else:
                target[neutral], target_genome[neutral] = seq, genome
        path = workdir / "proteomes" / f"g{gi:03d}.faa"
        path.parent.mkdir(exist_ok=True)
        path.write_text("".join(normalized))
        genome_records.append({"genome": genome, "role": role, "proteins": len(source),
                               "residues": sum(len(s) for _, s in source), "path": str(path),
                               "sha256": digest(path.read_bytes())})
    for filename, seqs in [("training.faa", train), ("queries.faa", target)]:
        (workdir / filename).write_text("".join(f">{p}\n{s}\n" for p, s in seqs.items()))
    dump(workdir / "protein-map.json", records)
    dump(workdir / "proteome-map.json", genome_records)
    return train, target, train_genome, target_genome, genome_records


def commands(policy, workdir):
    binary = policy["mmseqs_executable"]
    p = policy["search_parameters"]
    flags = ["-c", str(p["coverage"]), "--min-seq-id", str(p["minimum_sequence_identity_parameter"]),
             "-s", str(p["sensitivity"]), "--threads", str(p["threads"]), "-v", str(p["verbosity"])]
    tail = ["--threads", str(p["threads"]), "-v", str(p["verbosity"])]
    path = lambda name: str(workdir / name)
    return [
        [binary, "createdb", "stdin", path("train_db"), "-v", "3"],
        [binary, "cluster", path("train_db"), path("clusters"), path("tmp_cluster"), *flags],
        [binary, "createtsv", path("train_db"), path("train_db"), path("clusters"), path("clusters.tsv"), *tail],
        [binary, "search", path("train_db"), path("train_db"), path("train_hits"), path("tmp_self"), *flags],
        [binary, "createtsv", path("train_db"), path("train_db"), path("train_hits"), path("self.tsv"), *tail],
        [binary, "createdb", "stdin", path("query_db"), "-v", "3"],
        [binary, "search", path("query_db"), path("train_db"), path("query_hits"), path("tmp_query"), *flags],
        [binary, "createtsv", path("query_db"), path("train_db"), path("query_hits"), path("query.tsv"), *tail],
    ]


def freeze_vocabulary(policy, workdir, train, train_genome):
    families = read_families(workdir / "clusters.tsv", train)
    hits, _ = best_hits(workdir / "self.tsv", train, train, families)
    if set(hits) != set(train):
        raise ValueError("Training self-search did not cover every training protein")
    fmap, occurrence, _ = consolidate(policy["training_phages"], train_genome, hits)
    if not fmap:
        raise ValueError("No retained training features")
    dump(workdir / "feature-map.json", [{"family": f, "feature": fmap[f],
                                         "training_occurrence": list(occurrence[f])} for f in fmap])
    dump(workdir / "member-family-map.json", families)
    dump(workdir / "best-self.json", hits)
    names = ["training.faa", "clusters.tsv", "self.tsv", "feature-map.json", "member-family-map.json", "best-self.json"]
    freeze = {"completed_command_count": 5, "target_search_started": False,
              "training_genomes": policy["training_phages"],
              "artifacts": [{"path": str(workdir / n), "sha256": digest((workdir / n).read_bytes())} for n in names]}
    dump(workdir / "vocabulary-freeze.json", freeze)


def verify_freeze(workdir):
    freeze = json.loads((workdir / "vocabulary-freeze.json").read_text())
    if freeze["completed_command_count"] != 5 or freeze["target_search_started"]:
        raise ValueError("Invalid pre-target vocabulary boundary")
    for item in freeze["artifacts"]:
        checked(item["path"], item["sha256"])


def run_commands(cmds, workdir, timeout, after_training):
    ledger = []
    for index, command in enumerate(cmds):
        if index >= 5:
            verify_freeze(workdir)
        print(f"START step {index + 1}/{len(cmds)} {command[1]}", flush=True)
        log = workdir / f"step-{index + 1:02d}.log"
        usage = resource.getrusage(resource.RUSAGE_CHILDREN)
        start = time.monotonic()
        stdin_path = workdir / ("training.faa" if index == 0 else "queries.faa") if index in (0, 5) else None
        with log.open("wb") as stream:
            stdin = stdin_path.open("rb") if stdin_path else None
            child = subprocess.Popen(command, stdin=stdin, stdout=stream, stderr=subprocess.STDOUT, start_new_session=True)
            try:
                code = child.wait(timeout=timeout)
            except subprocess.TimeoutExpired:
                os.killpg(child.pid, signal.SIGKILL)
                child.wait()
                code = -signal.SIGKILL
            finally:
                if stdin:
                    stdin.close()
        end = resource.getrusage(resource.RUSAGE_CHILDREN)
        ledger.append({"command": command, "exit_code": code, "log": str(log),
                       "stdin": str(stdin_path) if stdin_path else None,
                       "stdin_sha256": digest(stdin_path.read_bytes()) if stdin_path else None,
                       "log_sha256": digest(log.read_bytes()), "wall_seconds": time.monotonic() - start,
                       "child_user_cpu_seconds": end.ru_utime - usage.ru_utime,
                       "child_system_cpu_seconds": end.ru_stime - usage.ru_stime,
                       "cumulative_max_child_rss_bytes": end.ru_maxrss * (1 if platform.system() == "Darwin" else 1024)})
        dump(workdir / "command-ledger.json", ledger)
        if code:
            raise RuntimeError(f"MMseqs failed or timed out at step {index + 1}: exit {code}")
        if index == 4:
            after_training()
    return ledger


def postprocess(policy, workdir, train, target, train_genome, target_genome, genomes):
    verify_freeze(workdir)
    families = read_families(workdir / "clusters.tsv", train)
    self_hits, self_summary = best_hits(workdir / "self.tsv", train, train, families)
    if set(self_hits) != set(train):
        raise ValueError("Training self-search did not cover every training protein")
    family_features, occurrence, counts = consolidate(policy["training_phages"], train_genome, self_hits)
    if not family_features:
        raise ValueError("No retained features")
    target_hits, target_summary = best_hits(workdir / "query.tsv", target, train, families)
    names = policy["training_phages"] + policy["target_phages"] + ["mapping_control"]
    features, vectors = project(names, train_genome | target_genome, self_hits | target_hits, family_features)
    mapping = [{"family": f, "feature": family_features[f], "training_occurrence": list(occurrence[f])}
               for f in family_features]
    if (mapping != json.loads((workdir / "feature-map.json").read_text())
            or self_hits != json.loads((workdir / "best-self.json").read_text())
            or families != json.loads((workdir / "member-family-map.json").read_text())):
        raise ValueError("Recomputed training vocabulary differs from frozen pre-target maps")
    dump(workdir / "best-query.json", target_hits)
    with (workdir / "feature-matrix.tsv").open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(["phage", *features])
        for name in names:
            writer.writerow([name, *vectors[name]])
    all_hits = self_hits | target_hits
    protein_genome = train_genome | target_genome
    assigned = Counter(protein_genome[p] for p in all_hits)
    retained = Counter(protein_genome[p] for p, h in all_hits.items() if h["family"] in family_features)
    rows = []
    for g in genomes:
        name, total = g["genome"], g["proteins"]
        rows.append(g | {"assigned_proteins": assigned[name], "assigned_protein_fraction": assigned[name] / total,
                         "proteins_mapped_to_retained_features": retained[name],
                         "retained_feature_count": sum(vectors[name]), "zero_vector": not any(vectors[name])})
    control_name = policy["training_phages"][0]
    return {"training_genomes": len(policy["training_phages"]), "training_proteins": len(train),
            "clustered_families": len(set(families.values())), "complete_member_coverage": True,
            **counts, "self_search": self_summary, "query_search": target_summary,
            "proteomes": rows, "target_vectors_equal": vectors["T7"] == vectors["T3"],
            "mapping_control": {"source_training_phage": control_name,
                                "exact_feature_vector_match": vectors[control_name] == vectors["mapping_control"],
                                "training_control_nonzero": any(vectors[control_name]),
                                "clean_nonzero_mapping_control": bool(any(vectors[control_name]) and vectors[control_name] == vectors["mapping_control"]),
                                "interpretation": "Mapping control only; not specificity validation"}}


def execute(policy_path, workdir, replay=False):
    helper_sha = digest(Path(__file__).read_bytes())
    policy_bytes = policy_path.read_bytes()
    policy = json.loads(policy_bytes)
    aa, verified = prepare(policy)  # Authenticate every actual annotation FASTA before commands.
    if replay:
        prior = json.loads((workdir / "run-inputs.json").read_text())
        if prior["policy_sha256"] != digest(policy_bytes):
            raise ValueError("Replay policy differs")
        if prior["verified_annotation_files"] != verified:
            raise ValueError("Replay source annotation inventory differs")
        expected_commands = commands(policy, workdir)
        if (prior["commands"] != expected_commands or prior["commands_sha256"] !=
                digest(json.dumps(expected_commands, separators=(",", ":")).encode())):
            raise ValueError("Replay command plan differs")
    elif workdir.exists():
        raise ValueError("Preserve existing raw workdir; use --replay for completed raw output")
    else:
        workdir.mkdir(parents=True)
    # Replay verifies normalized inputs against the pre-command hashes, not regenerated outputs.
    if replay:
        for item in prior["normalized_artifacts"]:
            checked(item["path"], item["sha256"])
        train = {n: s for n, s in proteins((workdir / "training.faa").read_bytes())}
        target = {n: s for n, s in proteins((workdir / "queries.faa").read_bytes())}
        pmap = json.loads((workdir / "protein-map.json").read_text())
        selections = [(g, "train", aa[g]) for g in policy["training_phages"]]
        selections += [(g, "target", aa[g]) for g in policy["target_phages"]]
        selections += [("mapping_control", "control", aa[policy["training_phages"][0]])]
        expected_map, expected_train, expected_target = [], {}, {}
        for gi, (genome, role, source) in enumerate(selections):
            for pi, (original, sequence) in enumerate(source):
                neutral = f"g{gi:03d}p{pi:06d}"
                expected_map.append({"protein": neutral, "genome": genome, "role": role,
                                     "source_protein": original, "sequence_sha256": digest(sequence.encode()),
                                     "length": len(sequence)})
                (expected_train if role == "train" else expected_target)[neutral] = sequence
        if pmap != expected_map or train != expected_train or target != expected_target:
            raise ValueError("Replay normalized proteins differ from authenticated source proteins")
        train_genome = {p["protein"]: p["genome"] for p in pmap if p["role"] == "train"}
        target_genome = {p["protein"]: p["genome"] for p in pmap if p["role"] != "train"}
        genomes = json.loads((workdir / "proteome-map.json").read_text())
        ledger = json.loads((workdir / "command-ledger.json").read_text())
        if len(ledger) != 8 or any(c["exit_code"] for c in ledger):
            raise ValueError("Replay requires all eight completed commands")
        if [c["command"] for c in ledger] != expected_commands:
            raise ValueError("Completed command ledger differs from frozen command plan")
        for c in ledger:
            checked(c["log"], c["log_sha256"])
        completed = prior.get("completed_artifacts", [])
        if [c["path"] for c in completed] != [str(workdir / n) for n in ["clusters.tsv", "self.tsv", "query.tsv"]]:
            raise ValueError("Missing/invalid completed external artifact inventory")
        for item in completed:
            checked(item["path"], item["sha256"])
    else:
        train, target, train_genome, target_genome, genomes = normalized_inputs(policy, aa, workdir)
        cmds = commands(policy, workdir)
        input_paths = [workdir / n for n in ["training.faa", "queries.faa", "protein-map.json", "proteome-map.json"]]
        input_paths += [Path(g["path"]) for g in genomes]
        prior = {"policy_sha256": digest(policy_bytes), "commands": cmds,
                 "commands_sha256": digest(json.dumps(cmds, separators=(",", ":")).encode()),
                 "verified_annotation_files": verified,
                 "normalized_artifacts": [{"path": str(p), "sha256": digest(p.read_bytes())} for p in input_paths]}
        dump(workdir / "run-inputs.json", prior)
        ledger = run_commands(cmds, workdir, policy["limits"]["timeout_seconds_per_command"],
                              lambda: freeze_vocabulary(policy, workdir, train, train_genome))
        # Freeze completed external output hashes before any Python postprocessing.
        prior["completed_artifacts"] = [{"path": str(workdir / n), "sha256": digest((workdir / n).read_bytes())}
                                        for n in ["clusters.tsv", "self.tsv", "query.tsv"]]
        dump(workdir / "run-inputs.json", prior)
    report = {"schema_version": 1, "purpose": policy["purpose"], "status": "complete",
              "policy_path": str(policy_path), "policy_sha256": digest(policy_bytes),
              "script_sha256": helper_sha, "raw_workdir": str(workdir),
              "software_deviation": policy["software_deviation"],
              "mmseqs_executable": policy["mmseqs_executable"], "mmseqs_sha256": policy["mmseqs_sha256"],
              "commands": ledger, "commands_sha256": prior["commands_sha256"],
              "phenotype_tables_loaded": False, "interaction_models_fitted": False,
              "interaction_scores_generated": False, "training_vocabulary_fitted": True,
              "train_only_vocabulary_and_occurrence": True,
              "normalization": policy["protein_normalization"],
              "createdb_input": "Sorted deterministic concatenated protein FASTA supplied as stdin, matching upstream createdb stdin semantics",
              "vocabulary_frozen_after_command": 5,
              "unoverridden_defaults": policy["unoverridden_defaults"],
              "effective_parameter_logs": [c["log"] for c in ledger],
              "resources": {"platform": platform.platform(), "architecture": platform.machine(),
                            "rss_scope": "Cumulative maximum RSS of serial direct child processes, not aggregate workflow memory"},
              "verified_annotation_files": verified,
              "interpretation": policy["interpretation"], "full_training_ready": False}
    report.update(postprocess(policy, workdir, train, target, train_genome, target_genome, genomes))
    if digest(Path(__file__).read_bytes()) != helper_sha:
        raise ValueError("Helper changed during execution; preserve raw output and replay with a stable helper")
    if replay and (workdir / "initial-report.json").exists():
        initial_bytes = (workdir / "initial-report.json").read_bytes()
        initial_report = json.loads(initial_bytes)
        checked(workdir / "executed-helper.py", initial_report["script_sha256"])
        if initial_report["policy_sha256"] != digest(policy_bytes):
            raise ValueError("Retained execution report refers to a different policy")
        report["replay_provenance"] = {
            "external_commands_repeated": False,
            "original_execution_helper": str(workdir / "executed-helper.py"),
            "original_execution_helper_sha256": initial_report["script_sha256"],
            "original_report": str(workdir / "initial-report.json"),
            "original_report_sha256": digest(initial_bytes),
            "refinement": "Reject internal/double-terminal stops; strengthen source/command/raw replay authentication; no parameter or numeric retuning"}
    artifact_names = ["run-inputs.json", "command-ledger.json", "protein-map.json", "proteome-map.json", "vocabulary-freeze.json",
                      "training.faa", "queries.faa", "clusters.tsv", "self.tsv", "query.tsv",
                      "feature-map.json", "member-family-map.json", "best-self.json", "best-query.json", "feature-matrix.tsv"]
    report["artifacts"] = [{"path": str(workdir / n), "sha256": digest((workdir / n).read_bytes())}
                           for n in artifact_names]
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--policy", type=Path, default=Path("software/tipcraft/benchmarks/retro-t7-t3/genophi-phage-feature-pilot.json"))
    parser.add_argument("--workdir", type=Path, default=Path("data/raw/genophi-phage-features/pilot-v1"))
    parser.add_argument("--output", type=Path, default=Path("data/processed/genophi-phage-feature-pilot.json"))
    parser.add_argument("--matrix-output", type=Path, default=Path("data/processed/genophi-phage-feature-matrix.tsv"))
    parser.add_argument("--replay", action="store_true")
    args = parser.parse_args()
    if args.output.exists() or args.matrix_output.exists():
        parser.error("Preserve existing processed outputs; choose new paths")
    # Processed output paths must never overwrite inputs/raw artifacts.
    forbidden = [args.policy.resolve(), args.workdir.resolve()]
    policy = json.loads(args.policy.read_text())
    forbidden += [Path(x["path"]).resolve() for x in policy["inputs"]]
    forbidden += [Path(policy["annotation_workdir"]).resolve(), Path(policy["mmseqs_executable"]).resolve()]
    if args.output.resolve() == args.matrix_output.resolve():
        parser.error("Report and matrix paths must differ")
    for out in (args.output, args.matrix_output):
        if any(out.resolve() == p or p in out.resolve().parents for p in forbidden):
            parser.error("Output overlaps protected input/raw path")
    try:
        report = execute(args.policy, args.workdir, args.replay)
    except (ValueError, OSError, RuntimeError) as exc:
        print(f"Pilot incomplete: {exc}", file=sys.stderr)
        return 1
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.matrix_output.parent.mkdir(parents=True, exist_ok=True)
    args.matrix_output.write_bytes((args.workdir / "feature-matrix.tsv").read_bytes())
    report["processed_matrix"] = {"path": str(args.matrix_output), "sha256": digest(args.matrix_output.read_bytes())}
    dump(args.output, report)
    print(json.dumps({"status": report["status"], "features": report["features_after_exact_occurrence_consolidation"],
                      "control": report["mapping_control"], "target_vectors_equal": report["target_vectors_equal"]}), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
