"""Train-only fixed-parameter fitting components for retrospective evaluation.

This module supplies no split search, hyperparameter search or production run
specification. Callers must freeze inputs, relatedness groups and evaluation
partitions separately. ML dependencies are imported only when fitting/scoring.
"""
from collections import Counter
from dataclasses import asdict, dataclass
import hashlib
import json
import math
from numbers import Real
from pathlib import Path
import platform
import sys


@dataclass(frozen=True)
class Sample:
    row_id: str
    host_id: str
    phage_id: str
    host_group: str
    phage_group: str
    publication: str


@dataclass(frozen=True)
class FitSettings:
    rfe_iterations: int
    final_iterations: int
    depth: int
    learning_rate: float
    l2_leaf_reg: float
    n_selected_features: int
    rfe_step: int
    seed: int
    threads: int
    feature_scope: str = "paired"
    weight_method: str = "none"


def _json_bytes(value):
    return (json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n").encode()


def _hash(data):
    return hashlib.sha256(data).hexdigest()


def _write_new(path, value):
    with path.open("xb") as stream:
        stream.write(_json_bytes(value))


def _text(value):
    return (isinstance(value, str) and bool(value) and value == value.strip()
            and not any(ord(c) < 32 or ord(c) == 127 for c in value))


def _samples(samples, require_nonempty=True):
    values = list(samples)
    if require_nonempty and not values:
        raise ValueError("Sample metadata must be nonempty")
    rows, pairs, hosts, phages = set(), set(), {}, {}
    for sample in values:
        if not isinstance(sample, Sample) or not all(_text(v) for v in asdict(sample).values()):
            raise ValueError("Each sample requires valid nonempty string identifiers/groups/publication")
        if sample.row_id in rows or (sample.host_id, sample.phage_id) in pairs:
            raise ValueError("Duplicate row identifier or host/phage pair")
        rows.add(sample.row_id)
        pairs.add((sample.host_id, sample.phage_id))
        for mapping, identifier, group in [(hosts, sample.host_id, sample.host_group),
                                           (phages, sample.phage_id, sample.phage_group)]:
            if identifier in mapping and mapping[identifier] != group:
                raise ValueError("One genome identifier maps to conflicting relatedness groups")
            mapping[identifier] = group
    return values


def validate_partition(train_samples, test_samples, regime="both_unseen", publication_holdout=False):
    """Validate explicit partitions; never generate or relax a split.

    cold_host blocks host IDs/groups; cold_phage blocks phage IDs/groups;
    both_unseen blocks both axes. Publication overlap is always disclosed.
    """
    if regime not in {"cold_host", "cold_phage", "both_unseen"}:
        raise ValueError("Unknown evaluation regime")
    if not isinstance(publication_holdout, bool):
        raise ValueError("publication_holdout must be boolean")
    train, test = _samples(train_samples), _samples(test_samples)
    # Also reject duplicate rows/pairs and inconsistent group assignments across partitions.
    _samples(train + test)
    overlaps = {}
    for name in ["host_id", "phage_id", "host_group", "phage_group", "publication"]:
        overlaps[name] = sorted({getattr(s, name) for s in train} & {getattr(s, name) for s in test})
    blocked = []
    if regime in {"cold_host", "both_unseen"}:
        blocked += ["host_id", "host_group"]
    if regime in {"cold_phage", "both_unseen"}:
        blocked += ["phage_id", "phage_group"]
    if publication_holdout:
        blocked.append("publication")
    violations = [name for name in blocked if overlaps[name]]
    if violations:
        raise ValueError("Partition overlap in protected axes: " + ", ".join(violations))
    return {"regime": regime, "publication_holdout": publication_holdout,
            "training_rows": len(train), "test_rows": len(test), "overlaps": overlaps,
            "publication_overlap_disclosed": overlaps["publication"], "random_fallback_used": False}


def _settings(settings):
    if not isinstance(settings, FitSettings):
        raise ValueError("settings must be FitSettings")
    for field in ["rfe_iterations", "final_iterations", "depth", "n_selected_features", "rfe_step", "threads"]:
        value = getattr(settings, field)
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            raise ValueError(f"{field} must be a positive integer")
    if settings.depth > 16:
        raise ValueError("CPU CatBoost depth must be at most 16")
    if isinstance(settings.seed, bool) or not isinstance(settings.seed, int) or not 0 <= settings.seed < 2**32:
        raise ValueError("seed must be an integer in [0, 2**32)")
    for name in ["learning_rate", "l2_leaf_reg"]:
        value = getattr(settings, name)
        if isinstance(value, bool) or not isinstance(value, Real) or not math.isfinite(value):
            raise ValueError(f"{name} must be finite numeric")
    if not 0 < settings.learning_rate <= 1 or settings.l2_leaf_reg < 0:
        raise ValueError("Require learning_rate in (0, 1] and nonnegative l2_leaf_reg")
    if settings.feature_scope not in {"paired", "host_only", "phage_only"}:
        raise ValueError("Unknown feature_scope")
    if settings.weight_method not in {"none", "per_phage_log10"}:
        raise ValueError("Unknown weight_method")


def _feature_names(feature_names):
    names = list(feature_names)
    reserved = set(Sample.__dataclass_fields__) | {"label", "labels", "target", "phenotype", "interaction"}
    if not names or any(not _text(n) or not n.startswith(("host:", "phage:"))
                        or not n.split(":", 1)[1] or n.split(":", 1)[1] in reserved for n in names):
        raise ValueError("Features require host:/phage: names; row metadata and labels are prohibited")
    if len(set(names)) != len(names):
        raise ValueError("Duplicate feature names")
    return names


def _eligible_indices(names, scope):
    if scope == "paired":
        return list(range(len(names)))
    prefix = {"host_only": "host:", "phage_only": "phage:"}[scope]
    return [i for i, name in enumerate(names) if name.startswith(prefix)]


def _matrix(X, width):
    if isinstance(X, (str, bytes, dict)):
        raise ValueError("Expected a two-dimensional binary matrix")
    rows = []
    try:
        for row in X:
            if isinstance(row, (str, bytes, dict)):
                raise ValueError("Expected binary matrix rows")
            values = list(row)
            if len(values) != width:
                raise ValueError("Feature matrix width differs from schema")
            if any(not isinstance(v, Real) or not math.isfinite(v) or v not in (0, 1) for v in values):
                raise ValueError("Matrix values must be finite numeric binary values")
            rows.append([int(v) for v in values])
    except TypeError as exc:
        raise ValueError("Expected a two-dimensional binary matrix") from exc
    if not rows:
        raise ValueError("Feature matrix must have rows")
    return rows


def _labels(y, count, both_classes=False):
    if isinstance(y, (str, bytes, dict)):
        raise ValueError("Expected a binary label vector")
    try:
        labels = list(y)
    except TypeError as exc:
        raise ValueError("Expected a binary label vector") from exc
    if len(labels) != count or any(not isinstance(v, Real) or not math.isfinite(v) or v not in (0, 1) for v in labels):
        raise ValueError("Label vector must match rows and contain finite binary values")
    labels = [int(v) for v in labels]
    if both_classes and set(labels) != {0, 1}:
        raise ValueError("Training requires both label classes")
    return labels


def _validate_genotype_features(matrix, samples, names):
    """A genome's representation cannot vary between interaction rows."""
    for prefix, identifier in [("host:", "host_id"), ("phage:", "phage_id")]:
        columns = [i for i, name in enumerate(names) if name.startswith(prefix)]
        observed = {}
        for row, sample in zip(matrix, samples):
            genome = getattr(sample, identifier)
            vector = tuple(row[i] for i in columns)
            if genome in observed and observed[genome] != vector:
                raise ValueError(f"Inconsistent {prefix} feature vector for repeated {identifier}")
            observed[genome] = vector


def training_weights(y_train, samples, method="none"):
    """Pure training-only weight rule; one-class vectors allowed here for audit.

    fit_bundle separately rejects a one-class training set. For each phage,
    positives receive max(1, 1+log10(nneg/(npos+1))); negatives receive 1.
    All-positive/all-negative phages receive unit weights.
    """
    sample_values = _samples(samples)
    labels = _labels(y_train, len(sample_values))
    if method not in {"none", "per_phage_log10"}:
        raise ValueError("Unknown weight method")
    counts = Counter(s.phage_id for s in sample_values)
    positives = Counter(s.phage_id for s, y in zip(sample_values, labels) if y)
    weights = []
    for sample, label in zip(sample_values, labels):
        npos = positives[sample.phage_id]
        nneg = counts[sample.phage_id] - npos
        weight = 1.0
        if method == "per_phage_log10" and label and npos and nneg:
            weight = max(1.0, 1.0 + math.log10(nneg / (npos + 1)))
        weights.append(weight)
    return weights


def _params(settings, iterations):
    return {"iterations": iterations, "depth": settings.depth,
            "learning_rate": settings.learning_rate, "l2_leaf_reg": settings.l2_leaf_reg,
            "random_seed": settings.seed, "thread_count": settings.threads,
            "loss_function": "Logloss", "task_type": "CPU", "verbose": False,
            "allow_writing_files": False, "use_best_model": False}


def fit_bundle(X_train, y_train, samples, feature_names, settings, output_dir, provenance=None):
    """Fit RFE and a fixed final estimator using only explicitly supplied rows.

    No evaluation data is accepted. This component does not certify provenance,
    partitions, preprocessing independence, or readiness for biological fitting.
    """
    _settings(settings)
    names = _feature_names(feature_names)
    matrix = _matrix(X_train, len(names))
    metadata = _samples(samples)
    if len(metadata) != len(matrix):
        raise ValueError("Sample metadata count differs from training rows")
    _validate_genotype_features(matrix, metadata, names)
    labels = _labels(y_train, len(matrix), both_classes=True)
    eligible = _eligible_indices(names, settings.feature_scope)
    if not eligible or settings.n_selected_features > len(eligible):
        raise ValueError("n_selected_features exceeds eligible feature count")
    weights = training_weights(labels, metadata, settings.weight_method)
    if provenance is not None and not isinstance(provenance, dict):
        raise ValueError("provenance must be a JSON object or None")
    provenance_bytes = _json_bytes(provenance or {})
    directory = Path(output_dir)
    if directory.exists():
        raise ValueError("Bundle output directory must be new")
    # Imports occur only after dependency-free validation, never in core CI.
    import numpy as np
    import sklearn
    from sklearn.feature_selection import RFE
    import catboost
    from catboost import CatBoostClassifier

    versions = {"python": sys.version.split()[0], "numpy": np.__version__,
                "scikit_learn": sklearn.__version__, "catboost": catboost.__version__}
    code_hash = _hash(Path(__file__).read_bytes())
    rfe_params, final_params = _params(settings, settings.rfe_iterations), _params(settings, settings.final_iterations)
    spec = {"schema_version": 1, "settings": asdict(settings), "feature_names": names,
            "eligible_global_indices": eligible, "training_rows": len(matrix),
            "input_hashes": {"matrix": _hash(_json_bytes(matrix)), "labels": _hash(_json_bytes(labels)),
                             "sample_metadata": _hash(_json_bytes([asdict(s) for s in metadata])),
                             "feature_schema": _hash(_json_bytes(names)), "weights": _hash(_json_bytes(weights))},
            "implementation_sha256": code_hash, "versions": versions,
            "rfe_parameters": rfe_params, "final_parameters": final_params,
            "provenance_sha256": _hash(provenance_bytes),
            "training_boundary": "Only supplied training rows/labels; no eval_set, early stopping or parameter search",
            "production_run_spec_complete": False,
            "caller_provenance_verified": False}
    directory.mkdir(parents=True, exist_ok=False)
    _write_new(directory / "fit-spec.json", spec)  # Freeze before either fit.
    with (directory / "provenance.json").open("xb") as stream:
        stream.write(provenance_bytes)
    spec_hash = _hash((directory / "fit-spec.json").read_bytes())
    with (directory / "fit-spec.sha256").open("x") as stream:
        stream.write(spec_hash + "\n")
    data, target, sample_weight = np.asarray(matrix, dtype=np.uint8), np.asarray(labels, dtype=np.uint8), np.asarray(weights, dtype=float)
    selector = RFE(CatBoostClassifier(**rfe_params), n_features_to_select=settings.n_selected_features,
                   step=settings.rfe_step)
    selector.fit(data[:, eligible], target, sample_weight=sample_weight)
    selected = [eligible[i] for i, keep in enumerate(selector.support_) if keep]
    if len(selected) != settings.n_selected_features:
        raise ValueError("RFE returned an unexpected feature count")
    selected_spec = {"global_indices": selected, "names": [names[i] for i in selected]}
    _write_new(directory / "selected-features.json", selected_spec)
    model = CatBoostClassifier(**final_params)
    model.fit(data[:, selected], target, sample_weight=sample_weight)
    model.save_model(str(directory / "model.cbm"), format="cbm")
    if _hash((directory / "fit-spec.json").read_bytes()) != spec_hash or _hash(Path(__file__).read_bytes()) != code_hash:
        raise ValueError("Fit specification or implementation changed during fitting")
    manifest = {"schema_version": 1, "status": "complete", "versions": versions,
                "artifacts": {name: {"sha256": _hash((directory / name).read_bytes())}
                              for name in ["fit-spec.json", "fit-spec.sha256", "provenance.json", "selected-features.json", "model.cbm"]},
                "selected_global_indices": selected, "selected_feature_names": selected_spec["names"],
                "rfe_effective_parameters": selector.estimator_.get_all_params(),
                "final_effective_parameters": model.get_all_params(),
                "feature_scope": settings.feature_scope, "weight_method": settings.weight_method,
                "weight_summary": {"minimum": min(weights), "maximum": max(weights), "rows": len(weights)},
                "execution_context": {"platform": platform.platform(), "architecture": platform.machine()},
                "score_type": "Uncalibrated CatBoost raw margin for class 1; not infection probability",
                "partition_validation": "Not performed here; caller must validate explicit partitions separately",
                "production_run_spec_complete": False, "held_out_data_arguments_accepted": False,
                "training_membership_verified": False}
    _write_new(directory / "manifest.json", manifest)
    return manifest


def score_bundle(X, feature_names, bundle_dir, *, expected_manifest_sha256):
    """Load authenticated CBM/schema artifacts and return finite raw margins.

    Accepts no phenotype vector, sample metadata or calibration outcomes.
    Callers freeze the manifest hash externally before scoring. This hash is
    checked before interpreting the manifest or any bundle schema/model.
    """
    directory = Path(bundle_dir)
    if (not isinstance(expected_manifest_sha256, str) or len(expected_manifest_sha256) != 64
            or any(c not in "0123456789abcdef" for c in expected_manifest_sha256)):
        raise ValueError("Expected a frozen lowercase SHA256 manifest digest")
    manifest_bytes = (directory / "manifest.json").read_bytes()
    if _hash(manifest_bytes) != expected_manifest_sha256:
        raise ValueError("Bundle manifest checksum differs from externally frozen digest")
    manifest = json.loads(manifest_bytes)
    names = _feature_names(feature_names)
    matrix = _matrix(X, len(names))
    if manifest.get("schema_version") != 1 or manifest.get("status") != "complete":
        raise ValueError("Incomplete or unsupported bundle")
    expected_files = {"fit-spec.json", "fit-spec.sha256", "provenance.json", "selected-features.json", "model.cbm"}
    if set(manifest.get("artifacts", {})) != expected_files:
        raise ValueError("Bundle artifact inventory differs")
    for filename in expected_files:
        if _hash((directory / filename).read_bytes()) != manifest["artifacts"][filename]["sha256"]:
            raise ValueError(f"Bundle checksum differs: {filename}")
    spec = json.loads((directory / "fit-spec.json").read_text())
    if (directory / "fit-spec.sha256").read_text() != _hash((directory / "fit-spec.json").read_bytes()) + "\n":
        raise ValueError("Frozen pre-fit specification checksum differs")
    selected = json.loads((directory / "selected-features.json").read_text())
    if spec.get("schema_version") != 1 or names != spec["feature_names"]:
        raise ValueError("Exact ordered feature schema required")
    settings = FitSettings(**spec["settings"])
    _settings(settings)
    indices = selected["global_indices"]
    eligible = _eligible_indices(names, settings.feature_scope)
    if (spec["eligible_global_indices"] != eligible or len(indices) != settings.n_selected_features
            or any(isinstance(i, bool) or not isinstance(i, int) or i not in eligible for i in indices)
            or len(set(indices)) != len(indices) or indices != sorted(indices)
            or selected["names"] != [names[i] for i in indices]
            or manifest["selected_global_indices"] != indices
            or manifest["selected_feature_names"] != selected["names"]
            or manifest["feature_scope"] != settings.feature_scope):
        raise ValueError("Selected feature map differs from fit specification")
    if _hash((directory / "provenance.json").read_bytes()) != spec["provenance_sha256"]:
        raise ValueError("Provenance differs from frozen fit specification")
    import numpy as np
    from catboost import CatBoostClassifier

    model = CatBoostClassifier()
    model.load_model(str(directory / "model.cbm"), format="cbm")
    values = model.predict(np.asarray(matrix, dtype=np.uint8)[:, indices], prediction_type="RawFormulaVal")
    scores = [float(v) for v in values]
    if len(scores) != len(matrix) or any(not math.isfinite(v) for v in scores):
        raise ValueError("Model returned nonfinite or missing scores")
    return scores
