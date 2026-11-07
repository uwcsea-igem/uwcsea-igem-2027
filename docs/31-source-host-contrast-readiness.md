# Independent source-host contrast readiness

Updated 2026-11-07. **The available source test cohort does not identify phage-dependent host ordering beyond a fixed host preference.** Of 49 group-disjoint candidate testing phages, six have complete literal `Yes`/`No` annotations for both hosts. Two provide mixed-state comparisons; both favor BW25113. This constrains the ranking validation claim, independently of model complexity or feature resolution.

## Frozen, descriptive audit

The [policy](../software/tipcraft/benchmarks/receptor-resource/source-host-contrast-policy.json) was recorded before calculating these distributions. The [script](../scripts/audit-receptor-host-contrasts.py) consumes the authenticated author `Table_S1_updated.tsv`, the [unchanged source-role partition](27-phage-groups-and-source-split.md) and the verified genomic group map. It retains all five buckets and all 260 cohort members; susceptibility annotations cannot move a phage between partitions.

Only exact `Yes` and `No` strings form complete pairs for this descriptive count. Every other literal is retained separately, including `Not assayed` and `Resistant`. This bookkeeping rule is not a validation that the annotations are comparable binary plaque labels. Production hosts and receptor annotations are not used, and the original Ando case file is not loaded. The source TSV includes target-related phages: their annotations are counted in the unchanged excluded bucket, never admitted to training. This is an annotation-aware source audit, not a claim of blindness or complete phenotype isolation. No model, score, calibration or threshold selection runs.

The [executed receipt](../data/processed/receptor-source-host-contrasts.json) binds source files, policy and implementation and retains every category's IDs and genomic groups.

| Literal pair, BW25113 / BL21 | Candidate training phages (groups) | Candidate testing phages (groups) |
| --- | ---: | ---: |
| Yes / Yes | 23 (12) | 4 (1) |
| Yes / No | 17 (7) | 2 (2) |
| No / Yes | 0 (0) | 0 (0) |
| No / No | 3 (3) | 0 (0) |
| Incomplete/nonbinary pair | 41 | 43 |
| All records, no partition changes | 84 | 49 |

Group counts are category-specific and must not be summed as independent groups: a genomic group may contain multiple annotation patterns. All 41 incomplete training pairs and 43 incomplete testing pairs have literal `Not assayed` for BL21. BW25113 annotations are populated for all of these records. Missing BL21 states cannot be treated as `No`, and complete-pair selection must not be described as full-cohort phenotype coverage.

Only **6/49 testing phages (12.2%)** have complete literal pairs; they occupy three genomic groups. The four Yes/Yes testing phages are in one group. The two mixed-state phages, `KCP072` and `NpM`, occupy two separate groups. These are completeness and grouping counts, not validated test observations admitted to a model.

## Identifiability consequence

For a two-host within-phage concordance endpoint, only opposite literal states create a directional comparison. Both observed test comparisons ask for BW25113 above BL21. A single fixed BW25113-above-BL21 ordering can satisfy every such comparison without examining the phage. The candidate training set's 17 complete mixed-state pairs have that same orientation. Consequently, a nominally perfect ranking result on this complete-pair test subset could not demonstrate an advantage over that fixed ordering.

This is a logical property of the inspected annotations, **not an executed baseline score or a model performance estimate**. No favorable ordering was encoded as a target predictor. Absence of the reverse orientation in complete pairs does not establish that the unassayed phages lack it.

Yes/Yes and No/No cases can still matter for a properly defined binary or quantitative endpoint. They provide no positive-negative ordering contrast within those two-host queries. A separate externally trained classification diagnostic might be possible after endpoint/stock/source admission, but it would require its own pre-score contract, appropriate group-level validation and matched baselines. It cannot silently become evidence for the original three-host ranking task or unseen-host generalization.

The [native Ando panel has a complementary attribution limitation](30-collision-endpoint-correction.md): only T3 contributes mixed-label native comparisons. Neither panel alone establishes a broadly phage-dependent host-ranking mechanism. Their records cannot be pooled or recoded to manufacture opposing directions, stock equivalence or apparent independent coverage.

## Decision for the next scoring run

**HOLD for using this source subset alone to validate phage-dependent ranking beyond a fixed host preference.** Do not fit a complex paired model and claim its ranking ceiling demonstrates interaction learning. A suitable independent evaluation source needs enough authenticated, group-separated observations to test the intended contrast, or an explicitly different scientifically justified endpoint with corresponding evidence.

The representation itself is not rejected by this label audit. The [family-conditioned representation](29-family-host-capacity.md) retains the native panel's required host contrasts. However, model admission also requires source-stock correspondence and valid independent phenotypes. Under the original strict unseen-host-group exclusions, both available source host backgrounds are still removed, leaving no training host domain; the [represented-host alternative](25-reference-host-overlap.md) remains a separate regime, not an implicit exemption.

No new original Ando case is covered. Independent-method coverage remains **0/10**, WIsH's negative result remains, and four exact exchange-construct cases still abstain. No binding, assembly or functional infection capability is inferred.

## Reproduction and checks

```bash
uv run --project software/tipcraft python scripts/audit-receptor-host-contrasts.py --output /tmp/receptor-source-host-contrasts.json
```

Five tests cover nonbinary literals, missing states, single versus opposing directions, related-group counting, duplicate IDs and missing mappings. Source IDs are never changed to increase the number of complete or informative comparisons.
