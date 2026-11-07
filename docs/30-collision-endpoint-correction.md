# Correcting representation gates to match the endpoint

Updated 2026-11-07. **Different hosts do not need different feature vectors when the evaluation does not require different outcomes for them.** The previous family-capacity interpretation imposed a stronger requirement than the recorded native ranking endpoint. The measured vectors remain unchanged; the decision consequence is corrected here.

## Authoritative comparison

The frozen [Ando case table](../software/tipcraft/benchmarks/retro-t7-t3/cases.csv) records T7 plaque formation on BL21, BW25113 and MG1655; native T3 plaque formation on BL21 and no detected plaques on BW25113 or MG1655 under the source conditions. Those qualitative observations have not changed. This comparison uses the repository's already extracted source observations. A fresh PMC retrieval during this correction returned a browser challenge, so it is not described as a new independent reading of the original figure.

The [existing endpoint contract](../software/tipcraft/benchmarks/retro-t7-t3/next-method-contract.json) defines within-phage positive-negative concordance, with ties worth 0.5, averaging queries with both labels. It requires finite scores for all three hosts per covered phage. It does **not** require different representations or scores for every pair of host names.

| Frozen representation | Host collision | Positive-negative comparisons forced to tie |
| --- | --- | ---: |
| Uniform-call global 5-mer union | All three hosts | 2, both for T3 |
| Family-conditioned binary | BW25113 and MG1655 | 0 |

The [executed audit](../data/processed/host-collision-endpoint-audit.json) reads the existing labels explicitly and checks these counts from the versioned vectors. Its [script](../scripts/audit-host-collisions-against-endpoint.py) generates no scores, fits no model, chooses no feature or threshold and adds no benchmark coverage. This is a **post hoc, label-aware interpretation correction**, not a new blinded representation test.

The global-union encoder still forces all informative comparisons to tie, supporting its particular ranking NO-GO. The family encoder separates the positive BL21 from both T3-negative hosts. Its remaining collision does not erase a required contrast in these six observations. That observation establishes neither correct learned ordering nor a prediction signal.

## What is preserved and what is superseded

The original family [policy](../software/tipcraft/benchmarks/receptor-resource/family-host-capacity-policy.json), [execution receipt](../data/processed/receptor-family-host-capacity.json), [replay receipt](../data/processed/receptor-family-host-capacity-replay.json), threshold, source inputs and raw artifacts remain unchanged and hash-bound. In particular, `all_three_targets_distinguishable` remains false. No new threshold search or feature run occurred.

The original policy's sentence calling all-three distinction a “necessary capacity diagnostic” was overstrict **as an admission condition for this endpoint**. This document supersedes that interpretation; it does not silently edit the executed policy. The finding that BW25113/MG1655 cannot receive different scores from these binary host inputs remains mathematically valid. Calling that finding alone grounds to stop this native ranking evaluation was the error.

## Attribution remains the decisive limitation

Only T3 contributes native positive-negative ranking comparisons: T7 has three positive observations and no negatives in this panel. An independently trained host-only model could therefore produce the same target ordering as a paired model. A native ranking improvement alone cannot demonstrate phage-dependent specificity. The existing contract already requires a matched host-only ablation and preserves that requirement.

This is why a distinct-host-vector check cannot replace independent phenotype evaluation and attribution controls. Conversely, diagnostic use of these known labels cannot be repackaged as a predictor that simply maps the BL21-like feature to susceptibility. Labels from this panel remain excluded from training, calibration and method parameters.

The current decision is **HOLD for independent model fitting/scoring pending source and training admission**, with the family representation not ruled out by these endpoint contrasts. Source-stock correspondence, the two represented host backgrounds, target-related group exclusions, label semantics and a frozen external validation regime still require resolution. The old both-unseen-group contract has not been relaxed into a known-host experiment. Original independent-method coverage remains **0/10**; no assembly, binding or functional-infection capability is inferred.

Five tests cover same-label collisions, mixed-label collisions, uninformative all-positive queries, duplicate observations, missing vectors and a multi-query panel that really does require all host pairs to separate. The complete 270-test suite passes with Biopython 1.86 enabled. To regenerate the diagnostic into a fresh output:

```bash
uv run --project software/tipcraft python scripts/audit-host-collisions-against-endpoint.py --output /tmp/host-collision-endpoint-audit.json
```
