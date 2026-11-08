# Authenticated grouped source panel and feature launch

For the subsequent first-task control failure and diagnosis, see [An exact-copy control detects a training/query encoding mismatch](41-training-query-identity-control-failure.md). The frozen launch record below is retained; the batch has not completed successfully.

Updated 2026-11-08. **The completed genomic comparison supports a fixed secondary panel of 50 phages and 23 hosts, with 1,093 observed pairs in nine crossed group blocks.** All 18 host/phage feature tasks were planned and hash-frozen before execution. This establishes actual source membership and preprocessing inputs; no interaction model or native prediction has been executed.

## Accepted genomic calculation

The [completed grouping receipt](../data/processed/official-variant-genome-groups.json) and [similarity matrix](../data/processed/official-variant-genome-similarities.tsv) cover 266 genome containers, retaining all previous nodes and the three additional official variants as bridges. The fixed 69.99% component rule produces 52 components across the full reference graph. That is not the number of independent source observations or the number of retained source-panel groups.

| Comparison | Result |
| --- | ---: |
| Unordered nonself pairs | 35,245 |
| Ordered nonself directions | 70,490 |
| Nonself count/rounded-similarity discrepancies | 0 |
| Original R calculation errors | 0 |
| Maximum unrounded similarity difference | 5.684341886080802e-14 percentage points |
| Predeclared absolute tolerance | 1e-12 percentage points |
| Dataset-specific self-controls | 266 / 266 passed |

The generic helper still records `status: discrepancy`, because its strict full-length self criterion fails for **TP15**: both implementations return 87,528 identical bases for a length of 87,728 with 200 ambiguous bases. This generic failure remains in the receipt. The predeclared dataset-specific rule requires exact R/Python self-count agreement for the four identified ambiguous references and full-length agreement for unambiguous references. All 266 satisfy that rule; it was not introduced after this result. The other three ambiguous references pass the generic length check too. No universal equivalence to the upstream workflow is claimed.

The complete search artifact remains explicitly qualified as recovered with an unknown search exit code, as recorded in the [alignment checkpoint](39-complete-common-context-alignment.md). The new Python and R calculator commands both have observed successful exits. No search, similarity threshold or source-selection rule was changed for this result.

## Exact secondary membership

The [frozen source panel](../data/processed/secondary-reference-panel.json) binds official reference genotypes, observed labels, measured groups, annotation artifacts and complete train/test/embargo memberships. Labels are stored separately from Sample metadata and never determine fold assignment.

| Quantity | Full source roster | Eligible panel |
| --- | ---: | ---: |
| Phage references | 52 | 50 |
| Host references | 25 | 23 |
| Observed binary spot-clearance pairs | 1,240 | 1,093 |
| Unobserved pairs, retained as unknown | 60 | 57 among eligible genomes |
| Host genetic groups | — | 15 |
| Phage genetic groups | — | 20 |

Excluded host components contain **ECOR13 and ECOR25**; excluded phage components contain **Bas64 and Bas65**. The mask excludes 147 observed rows, each with its reason retained. No missing observation is turned into a negative. The 50-by-23 eligible Cartesian universe has 1,150 pairs: 1,093 observed and 57 unknown. Later scoring coverage must distinguish complete coverage of eligible observations from the fraction of the full 1,240-observation source roster.

All three official variants receive measured genotype groups rather than inheriting the old author-cache genotype assignment. Their public reference correspondence is not proof of identity to the physical assay stocks; that limitation remains for all references. Spot clearance is not automatically plaque formation or productive infection.

The fixed seed-2027 assignment uses three host-group folds crossed with three phage-group folds. Mixed quadrants are embargoed; every retained observed pair is tested once. All nine training partitions contain both labels. Informative test-phage counts vary across blocks, which the later metric must report rather than treat as equal independent biological evidence. No split search or fold-count fallback followed these diagnostics.

## Execution boundary

The [18-task execution plan](../data/processed/secondary-feature-execution-plan.json), SHA256 `2daf06d281cab07095d8b724c6614f340596601e825d00523c72c2c9a92a09a1`, fixes the exact panel hash, command arguments, output destinations and per-task plan-only receipts. Those receipts bind normalized input hashes, selected genomes and all eight MMseqs commands before the first biological feature command. Tasks run serially with the existing four-thread parameters and no automatic retry.

The first host-axis task is running at this checkpoint. The batch retains command logs and compares executed normalized-input hashes and commands against the frozen plans. Execution/control failures stop the batch; empty vocabularies remain explicit fitting-readiness failures without a representation fallback. The batch prepares features, not a learned host ranking. Complete feature receipts and the full fitting/evaluation manifest still precede the first biological interaction-model fit.

Original Ando new-method coverage remains **0/10**. WIsH's negative result and the four exact exchange-construct abstentions remain unchanged. No binding, recognition-head, assembly or design capability is inferred from these genomic and preprocessing results.
