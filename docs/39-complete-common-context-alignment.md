# Complete common-context alignment evidence

For the subsequent completed numerical comparison, admitted panel and feature launch, see [Authenticated grouped source panel and feature launch](40-secondary-panel-and-feature-launch.md).

Updated 2026-11-08. **The augmented reference corpus now has a complete, validated common-database nucleotide alignment artifact.** All 273 record queries, all declared alignment-row counts and all record self-controls pass. This removes the missing-search-artifact prerequisite for the original-R-checked grouping calculation; it does not yet establish accepted groups, source admission or prediction performance.

The [official accession audit](35-basel-phage-accession-correspondence.md) found three official reference genotypes that differ from the author cache. The search retains all 263 previous genome containers as bridges and adds those three official genotypes as distinct nodes. Multi-record containers remain separate DNA records, giving 266 containers and 273 queries. Every query searched the same complete database under the original fixed parameters. No different database contexts are mixed, and no reference DNA is concatenated, repaired or designed.

## Accepted artifact evidence

The [recovery receipt](../data/processed/official-variant-alignments-recovered.json) records:

| Check | Result |
| --- | ---: |
| Complete unique record queries | 273 / 273 |
| Database records verified | 273 |
| Passed record self-controls | 273 / 273 |
| Completed batch query counts | 55, 56, 46, 47, 43, 26; sum 273 |
| Declared HSP rows, verified against actual rows | 1,251,570 |
| Original commented BLAST artifact | 2,252,691,968 bytes |
| Comment-free data-row artifact | 2,252,598,720 bytes |
| Query order retained; coordinate offsets applied at this stage | Yes; no offsets applied |

HSP counts are alignment records, not independent phages, infection observations or biological replicates. The receipt's `numeric_hsps` binding identifies the comment-free 14-field data-row file; that file also retains the aligned sequence columns. It is not an already calculated genome-similarity matrix.

Raw artifact SHA256: `be58e56e707b1a321708b53393931584393203f7894d1cc205742809f6ef07b8`.

Comment-free artifact SHA256: `e35366ce5a5b6f34a30f0e4990c29e0f0b3a3630f5649724691af4342cb57275`.

The original search was resumed, not rerun. Its monitoring parent had been lost, so the search's final process exit code remains **unknown**. The database-construction exit code is observed zero. The explicit [recovery protocol](../software/tipcraft/benchmarks/receptor-resource/official-variant-recovery-policy.json) accepts complete artifact evidence without inventing search-exit success; the original execution policy remains preserved. Query completion includes all six terminal batches and the final EOF trailer, not an assumption that the last trailer alone should name all 273 queries.

## Grouping calculation in progress

The fixed [grouping consumer](../scripts/group-official-basel-variants.py) has started under the separately selected [recovered-artifact policy](../software/tipcraft/benchmarks/receptor-resource/official-variant-group-recovery-policy.json). It rechecks input bindings and completion, applies only recorded within-container coordinate offsets, and calculates 35,245 unordered container pairs. Its Python calculation has finished; the original R implementation is running at this checkpoint. No alternate similarity threshold, reduced corpus or source-outcome-selected subset has been substituted.

Acceptance still requires exact directional-count agreement, the frozen unrounded-similarity tolerance, rounded agreement and the declared dataset-specific self-control checks. Group edges remain at 69.99%, retaining all nodes as bridges and excluding target/direct-gp17-seed components according to the fixed rule. Until the consumer writes a passing completed report, no augmented grouping or exact eligible secondary-panel manifest is asserted.

No phenotype table or interaction model is used by this search/recovery/grouping stage. No new Ando score has been produced: new-method coverage remains **0/10**, WIsH's six-case negative result remains, and the four exact construct cases still abstain. These sequence-comparison artifacts supply no evidence of binding or tail-fibre assembly compatibility.
