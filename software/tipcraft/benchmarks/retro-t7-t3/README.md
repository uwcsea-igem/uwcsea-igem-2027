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

## Next evaluation, before any claim of feasibility

Choose and freeze an actual prediction method, inputs, thresholds and method/version provenance without consulting held-out outcomes. Do not use phage names, donor identity, paper labels or a manually copied lookup as predictive evidence. Group future splits by construct and publication; use additional independent positive and negative assembly evidence for calibration. Record structure versions, seeds and compute environment if a structural backend is introduced. Compare to the baseline on the same cases and report uncertainty. Do not tune and test on the same ten related observations.

`targets.py` only validates structure references; `design.py` and `order.py` explicitly raise NotImplementedError. `assembly.py` records supplied evidence and unknowns, not a prediction of compatibility. `rank.py` computes classification metrics and does not rank designs. No time-to-design or ordering claim is supported.
