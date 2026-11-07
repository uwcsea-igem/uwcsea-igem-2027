# Family-conditioned host representation capacity

Updated 2026-11-07. **The frozen family-conditioned binary representation separates BL21 from the two K-12 references, but BW25113 and MG1655 still collide.** It cannot distinguish all three candidate hosts as the sole host input to a fixed paired predictor. This measured boundary applies to this representation and these inputs; it is not a rejection of every family-aware method.

## Protocol and inputs

The [policy](../software/tipcraft/benchmarks/receptor-resource/family-host-capacity-policy.json) was recorded before the first biological execution. The [producer](../scripts/check-receptor-family-host-capacity.py) uses the five existing, hash-authenticated, uniform Prodigal 2.6.3 proteomes from the [previous control](26-host-encoder-capacity.md). Only author BL21 and Keio fit the vocabulary: 4,273 and 4,314 proteins. The original BL21, BW25113 and MG1655 references, with 4,234, 4,314 and 4,319 proteins, are projection-only. These author reference assets are independent of the unavailable Picard host-genome collection. No additional genome download occurred.

All proteins in this run pass the standard-20-residue check after removal of one terminal stop marker. Protein boundaries and source order are preserved; neutral identifiers identify computational records, not biological stock identity.

Pinned MMseqs2 18 clusters and self-searches the source proteins with coverage parameter 0.8, minimum sequence-identity parameter 0.4, sensitivity 7.5 and four threads. Logs preserve effective defaults. Best-hit selection follows the existing kernel's float32 identity/E-value ordering and raw-first ties. This is an explicit version adaptation from the paper-listed MMseqs2 15.6f452.

Within each self-assigned protein family, contiguous 5-mers occurring in at least two distinct source proteins are retained. Repetitions within one protein do not inflate this count. Members with identical source-genome occurrence tuples are merged. Thus a member is a **family/k-mer pair**, rather than a global k-mer. The training proteins, clusters, self-search, vocabulary and implementation are hash-frozen before the target query database is created.

Projection follows the [audited reached-family arithmetic](28-family-assignment-semantics.md), using the pinned upstream CLI threshold of **0.5**, without a threshold search. A member matches when its k-mer occurs in a query protein assigned to its family. Repeated query matches count once. Only vocabulary members from reached families enter each feature's denominator. Entirely unmatched genomes receive explicit zero rows as a declared adaptation. This implementation is not certified equivalent to the authors' complete workflow or exact paper invocation.

## Executed result

The [execution receipt](../data/processed/receptor-family-host-capacity.json) reports 4,249 source families and 1,167,642 retained family/k-mer members, consolidated into three features. Both source self-projections exactly reproduce their fitted occurrence vectors. The target search starts only after that control passes. All five genomes remain in the result, and none has an all-zero vector.

| Reference input | Shared-pattern feature f0 | BL21-pattern feature f1 | Keio-pattern feature f2 |
| --- | ---: | ---: | ---: |
| Author BL21, fitting control | 1 | 1 | 0 |
| Author Keio, fitting control | 1 | 0 | 1 |
| Original BL21, projection only | 1 | 1 | 0 |
| Original BW25113, projection only | 1 | 0 | 1 |
| Original MG1655, projection only | 1 | 0 | 1 |

For the BL21-pattern feature, original BL21 matches 608/608 reached members; BW25113 and MG1655 each match 8/95. For the Keio-pattern feature, BL21 matches 0/170, while both K-12 targets match 2,021/2,021. Shared-feature fractions differ slightly among references and all exceed the frozen threshold; these diagnostic fractions were not supplied as extra continuous predictor inputs or selected as a new scoring method.

All 8,587 training proteins have self-search hits. The query search covers 12,834 of 12,867 target proteins, leaving 33 without hits; those proteins and their genomes are not silently deleted from the input inventory. Neither search has a recorded cross-family best-hit tie. The complete eight-command run took about 77 seconds on this computer. No classifier, phenotype-based feature selection or infection scoring ran.

## Decision

Family conditioning changes the previous global-union result: BL21 no longer collides with the other two hosts. **BW25113/MG1655 discrimination remains NO-GO for this frozen binary host representation as the sole host input.** For a fixed phage, a predictor receiving only the same host vector cannot assign these two hosts different scores. More downstream model complexity cannot repair identical inputs without additional information.

This does not establish biological equivalence, exact assay-stock correspondence, infection accuracy, binding or assembly compatibility. The two source host backgrounds already overlap the original reference groups; the result cannot establish unseen-host-group performance. Additional host diversity, different representations or continuous match features would be new methods requiring their own rationale and pre-score protocol, rather than post hoc repairs declared successful here. No threshold or family-assignment retuning followed this result.

The appropriate next method decision is to resolve a representation and labelled host context that can support the original three-host task. The candidate [84/49 phage split](27-phage-groups-and-source-split.md) alone cannot supply the missing host axis. Training a classifier now would not resolve the demonstrated BW25113/MG1655 collision.

Original Ando independent-method coverage remains **0/10**. WIsH's measured failure and the four exact-construct abstentions remain unchanged. This capacity diagnostic supplies no new infection accuracy or assembly success claim.

## Verification and reproduction

Five focused tests check within-protein occurrence counting, family conditioning, self-projection, reached denominators, unmatched-genome retention, duplicate query matches and the inclusive threshold. The complete 265-test suite passes with Biopython 1.86 enabled.

An independent [replay helper](../scripts/replay-receptor-family-host-capacity.py) authenticates the bound inputs, reconstructs normalized FASTAs, rebuilds the training-only vocabulary and reproduces **both complete projection objects**, including every reported numerator, denominator and fraction, without rerunning MMseqs2. Its [receipt](../data/processed/receptor-family-host-capacity-replay.json) additionally hashes the raw query alignments and normalized query inputs. This checks deterministic bookkeeping; it does not independently validate the biological method.

With authenticated input caches and the pinned binary, reproduce into fresh destinations:

```bash
uv run --project software/tipcraft python scripts/check-receptor-family-host-capacity.py --workdir data/raw/receptor-resource/family-host-capacity-reproduction-1 --output /tmp/receptor-family-host-capacity.json
uv run --project software/tipcraft python scripts/replay-receptor-family-host-capacity.py --report /tmp/receptor-family-host-capacity.json --output /tmp/receptor-family-host-capacity-replay.json
```
