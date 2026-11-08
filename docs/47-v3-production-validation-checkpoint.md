# V3 rebuilt-vocabulary production checkpoint

Updated 2026-11-08. **The previously failing host partition now completes under the prospectively fixed symmetric projection procedure, with no mapping discrepancy and an authenticated rebuilt vocabulary.** The [initial host receipt](../data/processed/secondary-features-v3-host1-phage0-host.json) has nine observed zero command exits; exact executed commands and normalized-input hashes match the [frozen v3 plan](../data/processed/secondary-feature-execution-plan-v3.json). Every listed raw artifact hash is verified.

## Actual accepted feature tasks

| Task | Training genomes | Query genomes | Features | Same assigned control families | Both unassigned control proteins | Assignment discrepancies | Equal control vector |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| host1-phage0-host | 18 | 5 | 816 | 4,960 | 2 | 0 | Yes |
| host0-phage0-phage | 32 | 18 | 113 | 74 | 0 | 0 | Yes |

The host training projection assigns hits to 81,293 proteins and leaves two unassigned. Its vocabulary retains 5,686 families occurring in more than one training genome, consolidated into 816 features. The external query vocabulary boundary is verified at completed command 6, with query database creation not yet started. The vocabulary and best-training-projection hashes still agree with that boundary after the complete task. No genome in this task has a zero vector.

The [first completed phage task](../data/processed/secondary-features-v3-host0-phage0-phage.json) likewise passes all nine observed commands, exact preplan/normalized-input comparisons, raw artifact hashes and the pre-query vocabulary boundary. All 74 copied proteins have matching assigned families; there are no matched absences or zero-vector genomes in this task.

## What the repair changes

The two host control absences remain explicitly unassigned; they are not counted as successful family identifications. The prospective v3 mapping-outcome/vector contract passes because both sides use ordinary external queries and agree. The v2 report remains a failure under its original family-identification requirement.

For `host1-phage0-host`, the complete final feature-name list and every training/query genome vector are exactly equal to the preserved failed v2 report. The rebuilt v3 vocabulary still records the changed ordinary-query hit coverage. This observed equality is specific to this partition; it does not assert equivalence for other partitions or native queries. No biological model or infection score was used to select this repair.

## Batch and decision status

Two of 18 v3 tasks are accepted at this checkpoint. After initial validation passed, the remaining-task runner started the frozen second stage, with at most two tasks active at once. `host0-phage0-host` and `host0-phage1-host` are running. The initial-stage ledger snapshot, executed runner hash and child process IDs are retained in the raw batch ledger. No automatic retry occurs.

Complete feature admission is not inferred from these two tasks. All 18 actual accepted receipts and the full nine-block fitting/evaluation manifest still precede biological model fitting. No partial-block prediction metric is reported. Original Ando new-method coverage remains 0/10, the four exact exchange-construct cases abstain and the measured WIsH negative result remains preserved. Binding and tail-fiber assembly compatibility are unvalidated.
