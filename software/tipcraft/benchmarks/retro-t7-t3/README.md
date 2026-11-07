# Retrospective T7/T3 benchmark v0.1

Purpose: establish a reproducible evaluation contract before implementing any design model. Ten qualitative observations are manually curated from [Ando et al., 2015, Fig. 3 and the Swapping tail fibers section](https://pmc.ncbi.nlm.nih.gov/articles/PMC4785837/). Labels describe plaque detection in that study, not universal host susceptibility; no numeric EOP values have been inferred. `T7_T3_gp17` and `T3_T7_gp17` denote the full-gp17 exchanges in the publication. No sequences or experimental procedures are included.

The set contains four native positive observations, two native below-detection observations and four published-exchange observations. Each row has a stable ID, host, phage, label, source, locator and interpretation. Assembly labels remain unknown because these rows are not direct assembly assays. Related rows are not independent replicates. A second-person source/figure review remains pending.

## Run from repository root

```bash
uv run --project software/tipcraft tipcraft benchmark --cases software/tipcraft/benchmarks/retro-t7-t3/cases.csv --output data/processed/retro-t7-t3.json
```

The predeclared always-positive baseline gives accuracy 0.60 and balanced accuracy 0.50 on these ten cases. It does not use the labels to make predictions. The JSON reports 0 model coverage and unknown assembly status. This is a genuine baseline calculation on curated observations, not a successful structural benchmark.

To evaluate independently produced predictions:

```bash
uv run --project software/tipcraft tipcraft benchmark --cases software/tipcraft/benchmarks/retro-t7-t3/cases.csv --predictions /absolute/path/to/predictions.csv
```

Prediction CSV columns: `case_id,score,method`; score is a finite probability in [0,1], with positive classification at 0.5. Missing predictions are allowed but coverage is explicit; missing cases are never silently called failures. Duplicate and unknown IDs, empty data and invalid scores are rejected. Both class recalls are required for balanced accuracy. A partial-set result cannot pass the gate. Input hashes record exactly which files were used.

## Independent WIsH experiment

A real nonconstant sequence-scoring experiment is now complete. We selected [WIsH, Galiez et al. 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5870724/) and ran the upstream executable at pinned commit `10e3f11779e926fb6dbdb20197c248367a75887f`, order 8, alpha 16. Models are fitted only to the three host reference genomes; known infection labels and donor identities are not inputs. It produces a host ranking from mean conditional log likelihood, not an infection probability or assembly score.

It covers the six native observations and abstains on the four exchange cases because verified complete-genome references for those exact constructs were not established. The [follow-up source audit](../../../../docs/08-construct-source-audit.md) documents the bounded search and supplement access limits. It ranks BW25113 first for both T7 and T3. T3's documented positive ranks below both negatives. All six scores survive an actual rerun with anonymous filenames and neutral FASTA headers. The [December report](../../../../docs/07-december-gate-report.md) records the failed discrimination, matching-subset comparator, leakage checks and limitations. This does not validate a receptor-specific design pipeline.

### Reproduce model construction and scores

Requires Python 3.11+, uv, curl, cmake and a C++11 compiler. Downloads public genome references and the pinned GPL upstream source into ignored `data/raw/wish/`; no system package installation, Docker or GPU is needed. All commands run from the repository root:

```bash
uv run --project software/tipcraft python scripts/prepare-wish.py --manifest software/tipcraft/benchmarks/retro-t7-t3/wish-inputs.json --workdir data/raw/wish
uv run --project software/tipcraft tipcraft wish-score --manifest software/tipcraft/benchmarks/retro-t7-t3/wish-inputs.json --genomes data/raw/wish --executable data/raw/wish/pinned-build/WIsH --output-dir data/processed/wish/named
uv run --project software/tipcraft tipcraft wish-score --manifest software/tipcraft/benchmarks/retro-t7-t3/wish-inputs.json --genomes data/raw/wish --executable data/raw/wish/pinned-build/WIsH --output-dir data/processed/wish/anonymous --anonymize
uv run --project software/tipcraft tipcraft wish-report --manifest software/tipcraft/benchmarks/retro-t7-t3/wish-inputs.json --cases software/tipcraft/benchmarks/retro-t7-t3/cases.csv --scores data/processed/wish/named/scores.json --anonymous-scores data/processed/wish/anonymous/scores.json --output data/processed/retro-t7-t3-wish.json
```

The first two scoring runs do not read `cases.csv`. The fourth command evaluates frozen scores. Raw sequence hashes are checked before execution; the report records the exact score and case hashes. Recompiled binary hashes can differ across platforms even when output scores agree. `wish-report` can replay the checked-in outputs without downloading genomes or compiling WIsH. CI uses this offline path. Include `--manifest` to bind the score JSON files to the frozen manifest and adjacent `llikelihood.matrix` files; a mismatch is rejected. Without that optional argument the legacy command still runs, but artifact validation is marked `not_verified` and a passing metadata audit is labeled `metadata_only`. Failed controls remain failed.

`wish-inputs.json` records all five exact accessions, retrieval URLs and sequence hashes. Preserve alternative-run outputs in a different directory; do not overwrite the historical run if a reference or method changes. Source/figure curation awaits independent human review. This exploratory run selected its method after the labels were visible and is not a blinded validation.

### Next gate evidence

Any next method needs a documented connection to receptor-level specificity, independently sourced calibration examples, exact construct provenance and frozen held-out evaluation. Do not use a donor-identity lookup to count known swaps as successful prediction. Group splits by construct and publication and do not tune and test on the same ten related observations. Direct assembly assessment requires its own positive/negative evidence.

`targets.py` only validates structure references; `design.py` and `order.py` explicitly raise NotImplementedError. `assembly.py` records supplied evidence and unknowns. `rank.py` computes classification metrics and does not rank designs. No design-to-ordering or assembly-compatibility claim is supported.

### Next-method applicability and input audit

The [primary-source review](../../../../docs/09-next-method-feasibility.md) prioritizes GenoPHI conditionally and records why other strain/taxonomy methods cannot simply be transferred. The [actual pinned-input inventory](../../../../data/processed/genophi-input-audit.json) is reproducible with `scripts/audit-genophi-sources.py` using the commands in that review. It audits public training-source tables, not host predictions. The [next-method contract](next-method-contract.json) fixes target scope and evaluation rules while explicitly leaving mandatory training/model manifests unresolved. Status is **HOLD**, new-method coverage **0/10**, performance **not evaluated**. Do not combine that inventory with WIsH's 6/10 coverage or call the contract a completed preregistration.
