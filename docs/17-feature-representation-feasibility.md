# Genome features and recognition-head coverage

Updated 2026-11-07. **A protein-family vocabulary learned from 93 eligible source phages represents T7 and T3 differently, while neither exact reference gp17 receives a feature assignment.** This separates two questions for TipCraft: whole-genome host-prediction feasibility and representation of an interchangeable recognition head. The current evidence supports proceeding with the former's input work; it does not supply a head-comparison engine.

## Frozen inductive pilot

The [pilot specification](../software/tipcraft/benchmarks/retro-t7-t3/genophi-phage-feature-pilot.json) was recorded before target feature assignment. Membership comes from the 94-phage modelling roster, with `427_P2` excluded by the [sequence-relatedness audit](15-relatedness-audit.md). `T7_Portugal` is also excluded from the collection but was never among those 94 modelling IDs: the retained count is **93**, not 92. The [uniform-reannotation policy](../software/tipcraft/benchmarks/retro-t7-t3/genophi-preprocessing-policy.json) applies to source and target proteins alike. No interaction table or Ando outcome is a feature-pilot input.

Only retained source proteins enter the training database, clustering, self-search and occurrence consolidation. After these steps, the family membership, feature mapping, training occurrence vectors and self-search assignments are hashed in `vocabulary-freeze.json` **before creating or searching the target database**. The target queries contain T7, T3 and an independently named copy of the lexicographically first retained training proteome as a mapping control. Queries search all training proteins and project into the frozen vocabulary; they do not create new families or features.

The procedure follows the pinned [GenoPHI family/assignment implementation](https://github.com/Noonanav/GenoPHI/blob/f7e30abc1b0bbba1ae76b8f58a3363bbbc5cd54f/genophi/mmseqs2_clustering.py) and [novel-genome assignment](https://github.com/Noonanav/GenoPHI/blob/f7e30abc1b0bbba1ae76b8f58a3363bbbc5cd54f/genophi/workflows/assign_features_workflow.py), with declared adaptations: uniform gene calling, MMseqs2 **18-8cc5c** rather than paper-listed 15.6f452, exact occurrence-tuple equality rather than unchecked hash grouping, and explicit preservation of zero-vector genomes. Families must occur in more than one selected training genome before consolidation. No supervised feature selection or interaction model is fitted.

Explicit search/clustering parameters are `-c 0.8 --min-seq-id 0.4 -s 7.5 --threads 4`. Logs preserve effective defaults: coverage mode 0, sequence-identity mode 0, search E-value 0.001 and result cap 300; clustering mode 0 and initial alignment mode 3; search alignment mode 2. These modes differ by stage, so the identity parameter must not be described uniformly as exact counted amino-acid identity. [Official MMseqs2 explanation](https://github.com/soedinglab/MMseqs2/wiki#how-does-mmseqs2-compute-the-sequence-identity).

## Actual feature outputs

The [verified replay report](../data/processed/genophi-phage-feature-pilot-verified.json) binds every source protein file, executable, command, raw alignment, vocabulary and output. The [full numeric matrix](../data/processed/genophi-phage-feature-matrix-verified.tsv) contains all 93 training rows plus both targets and the control.

| Stage | Measured result |
| --- | --- |
| Source proteins | 9,677 across 93 phages |
| Clustered families | 2,123; every training protein occurs in exactly one family |
| Families retained after selected-genome occurrence filter | 1,527 |
| Features after exact occurrence consolidation | 297 |
| Training self-search | Every one of the 9,677 proteins has a returned hit |
| Ambiguous best hits across different families | None in this run |

| Query | Called proteins | Proteins with returned assignments | Active features | Zero vector? |
| --- | ---: | ---: | ---: | --- |
| T7 reference | 49 | 33 | 10 | No |
| T3 reference | 43 | 32 | 5 | No |
| Copied 409_P1 control | 48 | 48 | 17 | No |

T7 and T3 vectors differ. The control reproduces its training row exactly, with nonzero features. These results show that the declared assignment path can encode these references and recover a known mapping; they are not infection predictions, host discrimination or a generalization estimate. Distinct vectors do not establish that the differing features contain useful host-specific information.

## Which part of gp17 is represented?

A [supplementary assignment audit](../data/processed/genophi-gp17-feature-mapping.json) identifies the exact gp17 proteins by their previously verified sequence hashes. Each occurs once in its target proteome: T7 `NP_042005.1` (553 aa) and T3 `NP_523342.1` (558 aa). **Neither query has any returned raw alignment or retained feature under the frozen search settings.** The audit is explicitly **post hoc after feature outputs**; it changes no parameter, vocabulary or endpoint and launches no new search.

Whole-genome vector differences therefore cannot be presented as evidence that this vocabulary captures these recognition heads. In this representation an unassigned protein contributes no features. A future native-phage ranking signal could arise from other genomic features; it would not establish sensitivity to gp17 replacement. The observation also does not prove the absence of remote/domain homologs, lack of binding, incompatibility or inability to infect. A recognition-head method needs its own representation and evidence.

The [evaluation contract](../software/tipcraft/benchmarks/retro-t7-t3/next-method-contract.json) records a host-only ablation requirement before model fitting. A native T3 ranking win alone cannot establish phage-dependent specificity if a phage-agnostic host model produces the same ordering. Its existing concordance/WIsH screening condition remains unchanged; attribution is a separate question.

## Resources, integrity and reproduction

The eight serial commands took **50.50 seconds**; maximum direct-child RSS was **961,773,568 bytes** (about 0.96 GB). These measurements belong to the original **x86_64 launcher environment**. Subsequent ARM64/Python 3.11 replay rechecks files and recomputes feature bookkeeping without repeating MMseqs. The verified report separates external execution context from postprocessing context. This small phage-only run does not estimate memory for all-host clustering, RFE or CatBoost training; current annotation/feature checks do not require a GPU.

The [original report](../data/processed/genophi-phage-feature-pilot.json) and [intermediate replay](../data/processed/genophi-phage-feature-pilot-replay.json) remain unchanged. Their exact helpers are recoverable from Git commits `c348b26` and `60a8a27`. The current helper rejects internal/double-terminal stops and authenticates the original execution report on replay. All three feature matrices are byte-identical. This refinement changed input/provenance checks and resource-context reporting, not the searches or numeric result.

Run from the repository root with the earlier authenticated caches. The retained original report must exist at `data/raw/genophi-phage-features/pilot-v1/initial-report.json` (an exact copy of the versioned original report):

```bash
uv run --project software/tipcraft python scripts/check-genophi-phage-features.py --replay --execution-report-sha256 891e396e438136e3337b398beb0b39a40ff611c5493e79e13d7530d535df9959 --output /tmp/genophi-phage-features.json --matrix-output /tmp/genophi-phage-features.tsv
uv run --project software/tipcraft python scripts/check-genophi-gp17-feature-mapping.py --output /tmp/genophi-gp17-mapping.json
```

Use new output paths. A fresh search instead requires a new `--workdir` and omits `--replay` and its execution-report hash; the unchanged policy and executable are checked before starting. Do not treat a different platform's timings as the original measurement. Source sequences, protein strings, databases and alignments remain in ignored raw storage.

## Next decision evidence

Complete the full host transfer and ANI exclusions, then freeze retained labels/groups, feature preparation, model parameters, seeds and validation splits before fitting. Any external validation fold must rebuild its reference vocabulary and occurrence consolidation from that fold's training genomes only; this 93-phage vocabulary is not a reusable CV shortcut. New-method observation coverage remains **0/10** until actual host scores exist. WIsH's six scores and four exact-construct abstentions remain unchanged, and assembly remains a separate gate.

For a possible future secondary panel, the [BASEL source audit](../data/processed/basel-validation-source-audit.json) maps all **52 phage IDs** and **25 host IDs** in the published 1,240-pair matrix to listed accessions. The roster follows actual matrix IDs: metadata flags contain four extra phages, and the host mapping contains two extra strains. All 25 host IDs occur in the original training table. This is the method's already-published novel-phage/known-host spot-assay panel, reused retrospectively; accession versions, sequence identity and homology have not been checked. It has not replaced or been selected to tune the T7/T3 evaluation, and missing/unclear outcomes are not negatives.

```bash
uv run --project software/tipcraft python scripts/audit-basel-validation-sources.py --matrix data/raw/genophi-external-validation/BASEL_ECOR_interaction_matrix.csv --archive data/raw/genophi-audit/GenoPHI_processed_data.zip --training-table data/raw/genophi-audit/data/interaction_matrices/ecoli_interaction_matrix.csv --output /tmp/basel-validation-source-audit.json
```
