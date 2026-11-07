# Reference-host overlap and evaluation regime

Updated 2026-11-07. **The available author-host references represent the same two close genomic background groups as the three original target references.** A future model must disclose this overlap rather than describing every query as an unseen host. This finding concerns declared reference representations; experimental-stock identity remains unresolved.

## Executed reference comparison

The [policy](../software/tipcraft/benchmarks/receptor-resource/host-relatedness-policy.json) was fixed before comparing these five references. It reuses the existing project near-host rule: estimated ANI at least 99.5% and at least 90% alignment fraction in each genome, with 0.05 ANI-point and 0.5 AF-point boundary buffers. These are conservative grouping rules, not strain-identity definitions or estimator confidence intervals.

The [executed audit](../data/processed/receptor-host-relatedness-audit.json) uses the pinned skani 0.4.0 binary, the established estimator parameters and one thread. All **25 directed comparisons**, including five self-controls, were reported and passed the declared completeness/control requirements. Original source files and sequence hashes were authenticated before and after execution. No phenotype labels or model scores were loaded.

| Author reference → original target reference | Reported ANI (%) | AF in author / target (%) | Near-reference rule |
| --- | ---: | ---: | --- |
| BL21 → BL21 | 100.00 | 99.17 / 99.85 | Meets rule |
| BL21 → BW25113 | 98.97 | 92.89 / 91.43 | Outside rule |
| BL21 → MG1655 | 98.97 | 93.11 / 91.45 | Outside rule |
| Keio → BL21 | 98.97 | 91.33 / 93.58 | Outside rule |
| Keio → BW25113 | 100.00 | 99.88 / 100.00 | Meets rule |
| Keio → MG1655 | 100.00 | 100.00 / 100.00 | Meets rule |

All pair directions remain in the hash-bound raw matrix. Connected components retain every reference node, using an edge when either directed comparison meets the near or buffered-boundary rule. The resulting components are:

- Author BL21 reference; original BL21 reference.
- Author Keio reference; original BW25113 reference; original MG1655 reference.

Reported ANI of 100.00 does **not** mean byte-identical sequence, identical experimental stocks, or equal susceptibility. The [host-payload audit](../data/processed/receptor-host-payload-audit.json) already establishes that both author-host sequence hashes differ from the corresponding original references. ANI estimates summarize aligned genomic sequence at the tool's reported precision; structural differences, unaligned sequence and experimental context remain relevant.

## Implications for an independent prediction run

A strict unseen-host-group regime using this same near-host exclusion rule would remove **both available author host backgrounds** from training. This source alone therefore cannot provide a nonempty training domain under that regime. This is a concrete limitation of the available source/target combination, not a conclusion that every strain-prediction method fails.

A different, explicitly declared **unseen-phage / represented-host-background** experiment could remain investigable after stock/reference correspondence, label semantics, source admission and phage grouping are resolved. It must disclose host-group overlap and require disjoint phage groups. Genomic input for MG1655 would not make its background independently unseen under this rule. A change of evaluation regime requires its own pre-score contract; it cannot silently satisfy the original method-specific both-unseen contract.

For any paired model, comparison with matched-subset host-only or training-prevalence baselines is necessary to distinguish phage-dependent information from a fixed host ordering. Existing constant baselines and the original ten-case denominator remain available; no new baseline or target score was manufactured by this audit. Novel phage and novel host generalization are separate claims, and neither establishes tail-fibre assembly compatibility.

The current independent-method coverage of the ten original observations remains **0/10**. Four exact Ando exchange constructs remain source abstentions. [Whole-phage grouping](27-phage-groups-and-source-split.md) is now complete and original-R-checked. It supplies exclusion and validation boundaries, not an admitted prediction model.

## Reproduction

```bash
uv run --project software/tipcraft python scripts/audit-receptor-host-relatedness.py --workdir data/raw/receptor-resource/host-relatedness-reproduction-1 --output /tmp/receptor-host-relatedness-audit.json
```

Use fresh output paths and the authenticated source/target caches. The helper rejects duplicate reference paths/IDs and retains a HOLD result for missing comparisons instead of substituting zero or asserting a complete group assignment.
