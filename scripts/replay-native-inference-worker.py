"""Score an authenticated opaque six-row packet with three saved bundles.

This worker performs inference only. Its caller authenticates the completed
biological run, constructs the packet, preserves the external identity map and
compares returned margins. Filesystem isolation and end-to-end anonymity are
not established here.
"""
import argparse
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import sys

from tipcraft import learning

SCOPES = {"paired", "host_only", "phage_only"}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def score_packet(packet, scorer=learning.score_bundle):
    if (not isinstance(packet, dict)
            or set(packet) != {"schema_version", "row_ids", "feature_names", "matrix", "bundles"}
            or type(packet["schema_version"]) is not int or packet["schema_version"] != 1):
        raise ValueError("Exact inference-only packet schema required")
    expected_ids = [f"q{i:03d}" for i in range(6)]
    if packet["row_ids"] != expected_ids:
        raise ValueError("Exactly six ordered opaque row IDs required")
    names = learning._feature_names(packet["feature_names"])
    matrix = learning._matrix(packet["matrix"], len(names))
    if len(matrix) != 6 or not isinstance(packet["bundles"], dict) or set(packet["bundles"]) != SCOPES:
        raise ValueError("Six feature rows and exactly three scopes required")
    # Validate every binding before loading any model.
    for scope, binding in packet["bundles"].items():
        if (not isinstance(binding, dict) or set(binding) != {"path", "sha256"}
                or not isinstance(binding["path"], str)
                or Path(binding["path"]).name != "manifest.json"
                or sha(binding["path"]) != binding["sha256"]):
            raise ValueError("Frozen bundle manifest required: " + scope)
        manifest = json.loads(Path(binding["path"]).read_text())
        if manifest.get("feature_scope") != scope:
            raise ValueError("Bundle scope differs from packet scope")
    methods = {}
    for scope in sorted(SCOPES):
        binding = packet["bundles"][scope]
        values = scorer(matrix, names, Path(binding["path"]).parent,
                        expected_manifest_sha256=binding["sha256"])
        # The production scorer enforces finite scores and complete coverage.
        # Retain that contract even with an injected scorer in a technical test.
        import math
        if len(values) != 6 or any(type(v) not in (float, int) or not math.isfinite(v) for v in values):
            raise ValueError("Six finite margins required")
        methods[scope] = {"model_id": binding["sha256"], "scores": dict(zip(expected_ids, values))}
    for binding in packet["bundles"].values():
        if sha(binding["path"]) != binding["sha256"]:
            raise ValueError("Bundle manifest changed during replay")
    return {"schema_version": 1, "status": "completed_opaque_native_inference",
            "methods": methods, "rows_per_scope": 6,
            "model_refitted": False, "target_labels_supplied": False,
            "sequence_to_feature_replay_executed": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--packet", type=Path, required=True)
    parser.add_argument("--expected-packet-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("Use a fresh replay output")
    if sha(args.packet) != args.expected_packet_sha256:
        raise ValueError("Packet checksum differs before parsing")
    result = score_packet(json.loads(args.packet.read_text()))
    if sha(args.packet) != args.expected_packet_sha256:
        raise ValueError("Packet changed during replay")
    result["provenance"] = {
        "packet": {"path": str(args.packet), "sha256": args.expected_packet_sha256},
        "worker_sha256": sha(__file__), "learning_sha256": sha(learning.__file__),
        "python": sys.version.split()[0],
        "packages": {n: version(n) for n in ["catboost", "numpy", "pandas", "scipy", "scikit-learn"]}}
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps({"status": result["status"], "rows_per_scope": 6}))


if __name__ == "__main__":
    main()
