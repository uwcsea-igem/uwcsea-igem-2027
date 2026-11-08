# V3 symmetric production projection protocol

Updated 2026-11-08, before any v3 feature search or biological model fit. **V3 uses ordinary external projection for training proteins and query proteins alike.** The [complete ECOR12 diagnostic](45-complete-symmetric-projection-diagnostic.md) motivates this repair: it separates 4,960 matching assignments from two matched absences without changing a phenotype, candidate, split, model setting or score threshold.

## Versioned implementation

The [source feature producer](../scripts/prepare-secondary-features-v3.py), [final-source/native feature producer](../scripts/prepare-native-transfer-features-v3.py) and [native fitter](../scripts/fit-native-transfer-v3.py) are separate files. Historical implementations, v1/v2 policies, all failed receipts and the six accepted v2 feature tasks remain unchanged. No old output is relabelled as passing v3.

The new [command templates](../software/tipcraft/benchmarks/receptor-resource/mmseqs-search-template-policy-v3.json) contain nine commands. Training-only clustering is unchanged. An additional database is created from the same normalized training FASTA; its independent prefix is used as the training projection query against the original training target DB. Both training and external projection SEARCH commands retain the same coverage 0.8, minimum identity 0.4, sensitivity 7.5, computed-residue alignment mode 3 and four threads. Vocabulary construction uses ordinary training projection outcomes and freezes after command 6, before external query database creation at command 7. Legacy filenames `self.tsv` and `best-self.json` are retained as artifact names; reports explicitly identify ordinary external training projection.

The v3 control compares every copied protein and separately reports `same_assigned_family`, `both_unassigned` and assignment discrepancies. It requires complete unique comparisons, zero discrepancies and exact genome-vector equality. One-sided absence remains a failure. Matched absence is never counted as successful family identification. Missing training hits and zero-vector genomes remain explicit; an empty axis vocabulary is not ready for fitting. This is a prospective change to the control contract, not a retroactive waiver of the v2 requirement.

The [source feature policy](../software/tipcraft/benchmarks/receptor-resource/secondary-feature-policy-v3.json) and [native feature policy](../software/tipcraft/benchmarks/receptor-resource/native-transfer-feature-policy-v3.json) bind the new producers and preserved failure/diagnostic evidence. New versioned fitting/evaluation and source/native design policies bind the resulting chain. The learning, partition, protein-feature and ranking kernels are unchanged. Source/native scoring settings, RFE, group folds, masks, comparators, endpoints and screens exactly match their v2 parents. Method selection previously saw outcomes; the study remains exploratory retrospective work.

## Frozen execution

All 18 real plan-only invocations pass before a v3 feature command starts. The [complete feature execution plan](../data/processed/secondary-feature-execution-plan-v3.json), SHA256 `cd553f808b6eb6d8db293da2493add096f2e20db9864ed99014f0f074b3c534d`, fixes exact commands, normalized-input hashes and fresh destinations for every block/axis.

The previously failed `host1-phage0-host` partition is scheduled first to check a rebuilt production vocabulary. It must complete with the v3 mapping-outcome/vector control and a nonempty vocabulary before the remaining tasks start. This priority is based on the technical failure, not a phenotype result. There is no automatic retry or partial-block prediction metric.

The machine reports 14 logical CPUs and 38,654,705,664 bytes of memory. After that initial task passes, the plan permits at most two independent feature tasks at once, each keeping its serial nine-command/four-thread procedure. This changes batch scheduling rather than individual tool or scientific settings. If a task fails, no new task starts; already-running tasks finish with their receipts preserved. All 18 accepted actual reports and a separate nine-block fitting/evaluation manifest still precede the first biological model fit.

## Validation and current conclusion

The package suite passes 393 tests with Biopython 1.86. Tests cover the nine-command source/native orchestration, vocabulary freezing before external queries, matched missing assignments, one-sided absence, vector discrepancies and duplicate controls. Five locked-runtime CatBoost integration tests pass, including saved-model worker replays on synthetic data and the versioned native fitting component. No recent remote CI result is claimed.

These checks do not establish a completed v3 production vocabulary, host-prediction accuracy, reliable generalization, binding or tail-fiber assembly. New-method Ando coverage remains 0/10; the four exact exchange-construct cases still abstain and WIsH's measured negative result remains preserved. Biological fitting remains pending the full actual feature and experiment manifests.
