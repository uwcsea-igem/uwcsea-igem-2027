# Complete-copy replay and versioned projection correction

Updated 2026-11-08. **Actual-residue identity restores family-assignment consistency for all 4,704 proteins in the saved copied-proteome control.** This broader technical replay supports a prospectively recorded v2 projection procedure. It is not host-prediction accuracy, an independent new-genome validation, or a waiver of the [failed v1 feature task](41-training-query-identity-control-failure.md).

## Complete copied-proteome evidence

The [probe policy](../software/tipcraft/benchmarks/receptor-resource/mmseqs-control-panel-probe-policy.json), SHA256 `7c8e7b9feebe60801cb041963e4debf83328cdd92b79bdf544ef14fcbfd673a1`, was fixed before the [diagnostic helper](../scripts/probe-mmseqs-control-panel.py) ran. It selects every saved ECOR14 control member, verifies each actual normalized FAA copy against its training sequence, and preserves the full original sequence databases and saved prefilter candidates. Both contexts request alignment mode 3; no phenotype or native labels are supplied.

The [completed receipt](../data/processed/mmseqs-control-panel-probe.json) records:

| Check | Result |
| --- | ---: |
| Actual sequence copies verified | 4,704 / 4,704 |
| Training-context queries with hits | 4,704 |
| External-copy queries with hits | 4,704 |
| Retained alignment rows per context | 72,476 |
| Cross-family best-key tie queries per context | 4 |
| Family-assignment mismatches | 0 |
| Diagnostic command exits | Six observed zero exits |

This test reuses the original reference clusters and the existing copied proteome. It does not rebuild the entire training vocabulary or validate a separately chosen genome. Matching every copied member's family removes the recorded mismatch in this replay; a production vocabulary must still be rebuilt and pass its own family/vector control. The v1 failed receipt, its one-bit discrepancy and the 17 unstarted v1 tasks remain intact.

## Minimal v2 method change

The new [source design](../software/tipcraft/benchmarks/receptor-resource/secondary-model-design-policy-v2.json), [search template](../software/tipcraft/benchmarks/receptor-resource/mmseqs-search-template-policy-v2.json) and [feature policy](../software/tipcraft/benchmarks/receptor-resource/secondary-feature-policy-v2.json) add **`--alignment-mode 3` to the training and query SEARCH steps only**. Reference clustering retains its original settings. Consequently the unchanged clustering parameter 0.4 refers to its original estimated similarity criterion, while projection SEARCH now requests computed aligned-residue identity at 0.4. These are explicitly distinguished, not described as a single unchanged identity measurement or an exact published-method replication.

Gene calls, reference genotypes, source inclusion/exclusion masks, rare-family filtering, occurrence consolidation, float32 best-key parsing, raw-first ties, coverage/sensitivity settings, model feature caps, weights, tree parameters, seeds, group splits and screening thresholds are unchanged. No data or label was removed to improve a score, and no biological model score exists when this revision is recorded. Existing parameterized consumers use the new policy files without source-code changes. Every v2 policy binds its v1 parent and the failed/focal/complete-copy evidence; no v1 policy is overwritten.

The [v2 fitting](../software/tipcraft/benchmarks/receptor-resource/secondary-fitting-policy-v2.json) and [v2 source evaluation](../software/tipcraft/benchmarks/receptor-resource/secondary-evaluation-policy-v2.json) bind this revised feature procedure. The corresponding [native design](../software/tipcraft/benchmarks/receptor-resource/secondary-native-transfer-design-policy-v2.json), [native features](../software/tipcraft/benchmarks/receptor-resource/native-transfer-feature-policy-v2.json), [native fitting](../software/tipcraft/benchmarks/receptor-resource/native-transfer-fitting-policy-v2.json) and [native evaluation](../software/tipcraft/benchmarks/receptor-resource/native-transfer-evaluation-policy-v2.json) retain the same five original query references and same-subset WIsH comparison. Source performance does not select alternative native settings.

## Prospective execution and interpretation

The [18-task v2 plan](../data/processed/secondary-feature-execution-plan-v2.json), SHA256 `296a1426b86bfc5a2f9704bd1a7eda15a3e27e43c2711e1d63117be234c367fe`, fixes all nine unchanged source blocks, both axes, exact inputs/commands and fresh v2 destinations before any v2 feature command. Real plan-only invocations pass for all 18 tasks. Full production controls remain mandatory; no failed v1 output is relabelled or mixed into v2 admission. Biological interaction fitting still requires the complete feature and fitting/evaluation manifests.

A separately fixed [conditional-order diagnostic](../software/tipcraft/benchmarks/receptor-resource/conditional-ordering-diagnostic-policy.json) adds interpretation without changing the primary screen. A paired model can outperform a separately fitted host-only comparator while retaining a fixed host order, for example through additive host and phage scores or different selected host features. The [label-free diagnostic kernel](../software/tipcraft/src/tipcraft/conditional_ordering.py) therefore records selected host/phage feature counts and matched host-pair score differences, order states and reversals within each authenticated fitted model. It uses neither labels nor missing-cell imputation, applies no new significance/admission threshold, and never pools margins across models. This avoids treating a comparator improvement alone as proof of conditional host ordering or receptor biology.

The package suite passes **373 tests** with Biopython 1.86, including ten synthetic tests for that diagnostic. The previously recorded four synthetic CatBoost integration tests remain the model-execution evidence; no biological interaction model has been fitted. Complete v2 production-batch outcomes, source ranking performance and native scores remain pending. Original Ando new-method coverage is **0/10**, WIsH's measured failure is retained, and the four exact exchange-construct cases still abstain. No binding, assembly or sequence-design capability is established.

## First production checkpoint

The first complete v2 block now has passing [host-axis](../data/processed/secondary-features-v2-host0-phage0-host.json) and [phage-axis](../data/processed/secondary-features-v2-host0-phage0-phage.json) receipts. Both actual eight-command runs have observed zero exits; their executed normalized-input hashes and commands match the corresponding pre-execution plans.

| Axis | Training genomes | Query genomes | Consolidated features | Copied-control proteins | Family mismatches | Control vector equals source |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| Host | 15 | 8 | 691 | 4,704 | 0 | Yes |
| Phage | 32 | 18 | 113 | 74 | 0 | Yes |

All 69,656 host training proteins have self-search hits. Neither report has a zero-vector genome. These are completed feature/mapping results, not host-specificity scores or a requirement that every named host have a distinct vector. A repeated feature count does not imply that the v2 feature definitions or vectors equal the failed v1 result.

Two of 18 feature tasks are accepted at this checkpoint; the fixed batch has started `host0-phage1-host`. Global fitting readiness is not asserted from this first block. The other feature reports and the complete fitting/evaluation manifest still precede the first biological model fit. No partial-block prediction metric is reported and no v1 failed output is reused as passing v2 evidence.

## Second production checkpoint

The `host0-phage1` [host](../data/processed/secondary-features-v2-host0-phage1-host.json) and [phage](../data/processed/secondary-features-v2-host0-phage1-phage.json) tasks also complete successfully. All eight observed command exits per task are zero. Actual command arguments and normalized-input hashes match the pre-execution plans; the normalized inputs and every listed output artifact pass SHA256 verification.

| Axis | Training genomes | Query genomes | Consolidated features | Copied-control proteins | Control family/vector equality | Zero-vector genomes |
| --- | ---: | ---: | ---: | ---: | --- | --- |
| Host | 15 | 8 | 691 | 4,704 | Pass | None |
| Phage | 33 | 17 | 143 | 74 | Pass | Bas67, Bas68 |

The two phage query genomes with zero vectors remain in the fixed panel. Zero means that no retained training-derived feature is present under this procedure; it does not mean an empty proteome, absence of infection or failure of the copied-proteome control. Their feature representation carries no distinguishing phage signal in this block. They are neither removed nor replaced, and their eventual predictions and limitations must be reported. No model score or accuracy result exists at this checkpoint.

Four of 18 feature tasks are now accepted. The batch has started `host0-phage2-host`; complete feature readiness and the actual nine-block fitting/evaluation manifest still precede biological model fitting.
