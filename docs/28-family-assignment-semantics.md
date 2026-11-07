# Family-aware projection: admission requirements

Updated 2026-11-07. **The published family-aware assignment code is a different candidate from the rejected global-union host encoder. Its missing-family behavior must be specified before a biological capacity test.** A synthetic execution establishes the behavior below; it does not establish host prediction or rescue the previous encoder.

## Source and executed evidence

The probe authenticates the exact [GenoPHI k-mer assignment source](https://github.com/Noonanav/GenoPHI/blob/f7e30abc1b0bbba1ae76b8f58a3363bbbc5cd54f/genophi/workflows/kmer_assign_features_workflow.py) at commit `f7e30abc1b0bbba1ae76b8f58a3363bbbc5cd54f`, Git blob `8f65647d058a09167456122d5f91eecc7eddcb97`. Three unchanged function AST nodes execute with their required dependencies; other upstream module code, sequence searches and model fitting do not execute. The [policy](../software/tipcraft/benchmarks/receptor-resource/family-assignment-probe-policy.json) was recorded before this synthetic run. The [receipt](../data/processed/genophi-family-assignment-probe.json) binds source, probe, policy and dependency versions.

The upstream code joins query best-hit families to the stored family/k-mer mapping (lines 282–291), checks each mapped k-mer within its assigned query protein (lines 186–203), deduplicates matched cluster/k-mer entries and constructs a denominator from cluster IDs remaining after that join (lines 382–406). It then pivots those rows into genome features (lines 421–432).

One synthetic feature has two vocabulary members, `alpha/AAAAA` and `beta/CCCCC`. A deliberately diagnostic threshold of 0.75 separates one-of-one from one-of-two arithmetic. This threshold is not a biological parameter recommendation or an accuracy-tuned value.

| Synthetic query | Assigned-family match fraction | Actual upstream feature | Fixed two-member denominator illustration |
| --- | --- | ---: | ---: |
| Both families assigned; both k-mers present | 2/2 | 1 | 1 |
| Alpha assigned and matching; beta unassigned | 1/1 | 1 | 0 (1/2) |
| Alpha assigned but its k-mer absent | 0/1 | 0 | 0 |
| No family assigned | No denominator row | Genome omitted | 0 |

The fixed-denominator column illustrates a different rule; it is not a correction certified by the authors. In particular, the upstream function's matching fraction is conditional on reached family mappings. Missing families do not automatically become failed vocabulary members. A completely unassigned genome is omitted by this function. A later workflow might handle missing rows separately; this probe makes no claim about every caller or the paper's exact inference invocation.

## Consequence for the next candidate

The existing [global-union collision](26-host-encoder-capacity.md) remains a measured NO-GO for that particular encoder. The family-aware workflow cannot be declared equivalent to it: family assignment and the matching threshold can change projections. Conversely, no family-aware capacity or host-ranking result has yet been measured here.

A subsequent biological test must freeze these choices before generating its feature results:

1. Which training-only protein families, member k-mers and occurrence-pattern map constitute the vocabulary, including whether family-aware or global training features are used.
2. Which best-hit identity, coverage, ordering and cross-family tie rules assign query proteins. Existing MMseqs2 18 and the upstream paper's version must be distinguished.
3. Whether matching fractions use reached-family or complete-vocabulary denominators, with an explicit threshold and no target-label tuning.
4. How completely unmatched genomes, ambiguous residues, zero vectors and feature ordering are preserved and reported.
5. Whether self-projection reproduces training vectors before any target projection. A representation with failed controls cannot proceed to model fitting.

The first candidate should preserve upstream reached-family arithmetic and explicitly add all-zero rows for unassigned genomes as a declared adaptation, rather than silently switching denominators. A separate fixed-denominator method would need its own rationale and protocol. Neither should be fitted to Ando labels or selected because it produces the desired host order.

With only the two available source host backgrounds, genome-pattern consolidation still has at most three nonzero training patterns, even if k-mers have family context. That bound does not by itself prove a query collision. A new capacity test is necessary. Distinguishable vectors, if obtained, would establish representational capacity rather than infection accuracy; source-stock correspondence, endpoint labels and independent model evaluation remain separate gates. The [candidate 84/49 phage partition](27-phage-groups-and-source-split.md) does not supply those missing host labels or prove generalization to unseen host groups.

## Reproduction

Place the unchanged source file at `data/raw/receptor-resource/genophi-assignment-source-v1/kmer_assign_features_workflow.py`, obtained from the pinned source link above. The script rejects a source hash or Git blob mismatch. Then run:

```bash
uv run --project software/tipcraft --with biopython==1.86 --with pandas==2.3.3 --with tqdm==4.67.1 python scripts/probe-genophi-family-assignment.py --output /tmp/genophi-family-assignment-probe.json
```

This uses synthetic proteins only and needs none of the unavailable Picard host-genome download. No independent infection scores were produced; original Ando independent-method coverage remains **0/10**, with the WIsH failure and four exact-construct abstentions preserved.
