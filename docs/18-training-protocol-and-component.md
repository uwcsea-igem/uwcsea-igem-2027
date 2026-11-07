# Training protocol and fitting component

Updated 2026-11-07. **TipCraft now has an executed train-only RFE/CatBoost fitting component, tested on synthetic data, while the biological training run remains unstarted.** The component makes feature selection, weights, model settings and bundle authentication inspectable. Full genome membership and preprocessing independence still require the actual experiment's source audits and run specification.

## Why the upstream orchestration needs adaptation

The [pinned source-flow audit](../data/processed/genophi-training-flow-audit.json) authenticates four implementation files and records exact call locations. In the [protein-family workflow](https://github.com/Noonanav/GenoPHI/blob/f7e30abc1b0bbba1ae76b8f58a3363bbbc5cd54f/genophi/workflows/protein_family_workflow.py#L239), family construction precedes repeated feature selection. Each RFE fits its own training rows; selected-feature occurrence is then pooled, reduced tables retain all labelled rows, and [later modelling re-splits those tables](https://github.com/Noonanav/GenoPHI/blob/f7e30abc1b0bbba1ae76b8f58a3363bbbc5cd54f/genophi/select_feature_modeling.py#L510). A later held-out row can therefore have influenced earlier pooled selection. Globally constructed vocabulary also uses supplied held-out genome content.

The inspected splitting code can fall back to row-random partitions, including when clustering is disabled or grouped splitting is invalid. Candidate selection and early stopping use the supplied test partition as tuning data. These paths do not provide this project's untouched external evaluation. CLI and lower-level defaults also differ, and some weighting/seed settings do not propagate through every call.

These are findings about the inspected wrapper paths. The audit did not establish the manuscript's actual commands or demonstrate Ando target-label leakage. It does not revise the paper's reported results. The adaptation is to retain useful representation/fitting primitives while explicitly controlling partitions, preprocessing and evaluation for this experiment.

## Implemented fitting boundary

The [component](../software/tipcraft/src/tipcraft/learning.py) provides:

- Immutable sample metadata and explicit `cold_host`, `cold_phage` and `both_unseen` partition checks. Genome IDs and relatedness groups are checked separately; invalid partitions fail without a random fallback. Publication overlap is disclosed and can be prohibited when the protocol claims publication holdout.
- Finite binary `host:`/`phage:` feature matrices. Repeated genome IDs must have consistent corresponding feature vectors, preventing a pair-specific response from masquerading as a host-only genomic feature.
- Training-only RFE, explicit per-fit seed/settings and optional per-phage log10 weights computed from training labels. Positive weights follow `max(1, 1 + log10(n_negative / (n_positive + 1)))`; negatives and single-class phage groups receive unit weights. A globally single-class fit fails.
- Fixed iteration counts with no evaluation-set, early-stopping or parameter-search argument. Host-only fitting removes phage columns before RFE.
- A fit specification and its SHA-256 sidecar written before fitting, followed by selected-feature JSON, a CatBoost CBM file and an artifact manifest. Scoring requires the caller's previously frozen manifest hash before reading the bundle. It accepts feature values/schema but no phenotype vector or sample metadata.

The returned scores are **uncalibrated CatBoost raw margins**. They are suitable for a declared ranking experiment; they are not physical binding scores or calibrated infection probabilities. The core package imports ML libraries lazily, with a separate [locked learning environment](../software/tipcraft/experiments/learning/README.md).

The partition validator and fitting function are separate primitives. Calling `fit_bundle` does not authenticate the biological membership of `X_train`, prove that ANI groups are correct, or establish that feature construction excluded held-out genomes. Manifests explicitly record `training_membership_verified: false` and `production_run_spec_complete: false`. A full biological runner must check those owning sources and apply partition validation before fitting; this component does not replace that work.

## Actual engineering verification

The [self-check report](../data/processed/learning-engine-check.json) records an actual run on a fixed **4-host × 8-phage synthetic fixture**, with six binary features. Each synthetic genome has a consistent feature vector. No real sequence, published interaction table or Ando outcome is accepted by the runner.

Three models are fitted: paired features, a bijective renaming of the same sample metadata, and a host-only ablation. Their manifest hashes are saved in an external ledger before any scoring. The run verifies finite nonconstant scores, save/load consistency, renaming invariance, host-only invariance to phage-column perturbations, protected-group rejection, ordered schema checks, and rejection of modified manifest/model/specification/selection files. Score comparisons use the fixed tolerance 10⁻¹². These are software integrity checks, with no reported biological accuracy or coverage gain.

The runtime uses Python 3.12.5, CatBoost 1.2.7, scikit-learn 1.5.2, NumPy 1.26.4, pandas 2.2.2 and SciPy 1.14.1, with remaining dependencies pinned by `uv.lock`. The [synthetic settings](../software/tipcraft/experiments/learning/synthetic-check-spec.json) use small engineering-test iteration counts; they are **not selected production hyperparameters**. The CI configuration now includes this isolated integration check; local execution has passed, while remote CI awaits a published branch/run.

```bash
uv run --project software/tipcraft/experiments/learning --locked python -m unittest discover -s software/tipcraft/experiments/learning/tests -v
uv run --project software/tipcraft/experiments/learning --locked python scripts/check-learning-fit.py --workdir data/raw/learning-engine-check/reproduction-1 --output /tmp/learning-engine-check.json
```

## Required full-run protocol

After the complete host cache and target exclusions are available, freeze exact retained labels, genome IDs, relatedness groups, partitions, annotation, feature construction, feature selection, weights, model parameters and seeds. For every outer fold, fit the vocabulary, occurrence consolidation and any supervised selection only on its training partition. Any tuning or early stopping must stay within inner training folds; the present fixed-parameter component supplies no tuning stage. Freeze each resulting pipeline before evaluating the outer fold once.

The final T7/T3 model may then fit eligible external training data and project the targets through its frozen feature vocabulary. Save the actual model/reference bundle hashes before target inference and perform the required neutral-name rerun and evaluation-only label permutations. The host-only comparison remains part of attribution. The [next-method contract](../software/tipcraft/benchmarks/retro-t7-t3/next-method-contract.json) continues to govern these decisions; no component test changes its **0/10 new-method coverage**, four construct abstentions or separate assembly gate.

The [source-admission implementation](20-source-admission-and-partitions.md) now supplies the upstream roster/group checks and fixed crossed-block algorithm. Its actual missing-input HOLD prevents final roster/partition output until complete host authentication and target-related exclusions are available. Fold-specific vocabularies and the full model run specification remain subsequent steps.

The [reusable feature kernel](17-feature-representation-feasibility.md#reusable-training-only-feature-kernel) now verifies training-only occurrence handling and reproduces the existing phage pilot matrix exactly. This supplies the bookkeeping component for future fold-specific vocabularies; each fold still requires its own authenticated clustering/search execution before fitting.
