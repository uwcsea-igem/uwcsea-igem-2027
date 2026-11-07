"""Group published host distances; no model fitting or phenotype prediction."""
import argparse
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import re

SOURCE_COMMIT = "99656d8c020fcd1dbbc30a02ee10274a07655562"
MATRIX_BLOB = "cccab7b1e343f33af8d0eaa759c89e6480980af0"
MODEL_SHA256 = "c4795d786a1e13ac998e263dbf4aaf95fc0f94c52b842d02b6f7fc6e99ccfef2"
MODEL_COMMIT = "f7e30abc1b0bbba1ae76b8f58a3363bbbc5cd54f"
THRESHOLD = 1e-4
SOURCE_PATH = "data/genomics/bacteria/isolation_strains/panacota/tree/370+host_distance_matrix_orignames.tsv"


def unique_names(names: list[str]) -> dict[str, str]:
    index = {}
    for name in names:
        key = re.sub(r"[^a-z0-9]", "", name.casefold())
        if not key or key in index:
            raise ValueError("Duplicate, empty or colliding normalized host identifiers")
        index[key] = name
    return index


def group_hosts(matrix_text: str, training_hosts: list[str]) -> dict:
    rows = list(csv.reader(io.StringIO(matrix_text), delimiter="\t"))
    if not rows or len(rows[0]) < 2 or rows[0][0] != "bacteria":
        raise ValueError("Expected a labelled host-distance matrix")
    columns = rows[0][1:]
    source_index = unique_names(columns)
    row_names = [row[0] for row in rows[1:] if row]
    unique_names(row_names)
    if len(rows) != len(columns) + 1 or set(row_names) != set(columns):
        raise ValueError("Distance matrix must be square and complete")
    distances = {}
    for row in rows[1:]:
        if len(row) != len(columns) + 1:
            raise ValueError("Incomplete distance matrix row")
        values = [float(x) for x in row[1:]]
        if any(not math.isfinite(x) or x < 0 for x in values):
            raise ValueError("Distances must be finite and nonnegative")
        distances[row[0]] = dict(zip(columns, values))
    for a in columns:
        if distances[a][a] != 0:
            raise ValueError("Distance diagonal must be zero")
        for b in columns:
            if distances[a][b] != distances[b][a]:
                raise ValueError("Distance matrix must be symmetric")
    training_index = unique_names(training_hosts)
    if set(training_index) - set(source_index):
        raise ValueError("Unmapped training host identifiers")
    source_to_training = {source_index[key]: name for key, name in training_index.items()}
    parent = {name: name for name in columns}
    def find(name):
        while parent[name] != name:
            parent[name] = parent[parent[name]]
            name = parent[name]
        return name
    ordered = sorted(columns)
    edge_count = 0
    for i, a in enumerate(ordered):
        for b in ordered[i + 1:]:
            if distances[a][b] < THRESHOLD:
                parent[find(a)] = find(b)
                edge_count += 1
    members = {}
    for name in ordered:
        members.setdefault(find(name), []).append(name)
    groups, source_groups = {}, {}
    for names in members.values():
        digest = hashlib.sha256(json.dumps(names, ensure_ascii=True, separators=(",", ":")).encode()).hexdigest()
        group_id = "host-group-" + digest[:16]
        if group_id in groups:
            raise ValueError("Group identifier hash collision")
        groups[group_id] = {"source_members": names,
                            "training_members": sorted(source_to_training[n] for n in names if n in source_to_training)}
        source_groups.update({name: group_id for name in names})
    target_flags = {}
    for target in ("BL21", "BW25113", "MG1655"):
        name = source_index.get(re.sub(r"[^a-z0-9]", "", target.casefold()))
        target_flags[target] = {"matched_source_id": name, "group_id": source_groups.get(name),
                                "status": "provisional_name_match" if name else "not_in_source_roster_by_name"}
    return {"source_hosts": len(columns), "mapped_training_hosts": len(training_hosts),
            "threshold": THRESHOLD, "edge_rule": "strictly distance < 1e-4",
            "edges_below_threshold": edge_count, "group_count": len(groups),
            "groups": dict(sorted(groups.items())),
            "training_host_groups": {source_to_training[n]: source_groups[n] for n in sorted(source_to_training, key=source_to_training.get)},
            "named_target_groups": target_flags}


def run(matrix: Path, model_table: Path) -> dict:
    raw, model = matrix.read_bytes(), model_table.read_bytes()
    blob = hashlib.sha1(b"blob " + str(len(raw)).encode() + b"\0" + raw).hexdigest()
    if blob != MATRIX_BLOB:
        raise ValueError("Host-distance matrix differs from pinned source")
    if hashlib.sha256(model).hexdigest() != MODEL_SHA256:
        raise ValueError("GenoPHI table differs from pinned source")
    rows = list(csv.DictReader(io.StringIO(model.decode())))
    if not rows or set(rows[0]) != {"strain", "phage", "interaction"}:
        raise ValueError("Unexpected model-table schema")
    result = group_hosts(raw.decode(), sorted({row["strain"] for row in rows}))
    return {"schema_version": 1, "source_commit": SOURCE_COMMIT,
            "source_url": f"https://raw.githubusercontent.com/mdmparis/coli_phage_interactions_2023/{SOURCE_COMMIT}/{SOURCE_PATH}",
            "matrix_git_blob_sha1": blob, "matrix_sha256": hashlib.sha256(raw).hexdigest(),
            "model_table_sha256": MODEL_SHA256,
            "model_table_url": f"https://raw.githubusercontent.com/Noonanav/GenoPHI/{MODEL_COMMIT}/data/interaction_matrices/ecoli_interaction_matrix.csv",
            "threshold_reference": f"https://github.com/mdmparis/coli_phage_interactions_2023/blob/{SOURCE_COMMIT}/README.md",
            "method": "Connected components over every host in the published core-genome distance matrix",
            "name_mapping": "Unique case-insensitive alphanumeric identifiers; no phenotype matching",
            **result,
            "original_cv_folds_recovered": False, "training_or_scoring_executed": False,
            "phenotype_values_used_for_grouping": False,
            "interpretation": "Reconstructed related-host groups, not original CV folds. All source hosts retained as possible bridges. Named target flags do not establish stock identity or replace target-reference genome/ANI exclusions; no automatic fold assignment."}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ("matrix", "model-table", "output"):
        parser.add_argument("--" + flag, type=Path, required=True)
    args = parser.parse_args()
    if args.output.resolve() in {args.matrix.resolve(), args.model_table.resolve()}:
        parser.error("Output cannot overwrite inputs")
    report = run(args.matrix, args.model_table)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
