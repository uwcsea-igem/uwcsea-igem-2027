# An exact-copy control detects a training/query encoding mismatch

Updated 2026-11-08. **The first biological feature task stopped on an exact-copy mapping control before any interaction model was fitted.** The same 385-aa protein acquired different reported identity values in training-self and external-query contexts, changing one binary feature. A prespecified focal diagnostic reproduces the asymmetry and removes it when actual aligned-residue identity is requested. The original failed task remains failed; the focal diagnostic does not validate the full feature pipeline.

## Observed failure and preserved boundary

The [failed host feature receipt](../data/processed/secondary-features-host0-phage0-host.json) completed all eight MMseqs commands with observed zero exits. Its training vocabulary used 15 genomes and 69,656 predicted proteins, all with self-search hits: 8,833 clusters, 5,781 retained families and 691 consolidated occurrence features. These counts are representation/preparation results, not infection predictions.

The cloned ECOR14 control contains 4,704 proteins. **One family assignment differs, and the resulting genome vector differs at `pc_318` (index 244), from training 0 to copied-query 1.** Exact-copy agreement is therefore not satisfied. The [batch failure receipt](../data/processed/secondary-feature-batch-failure.json) records one started task, zero validated completed tasks and 17 unstarted tasks. The command exit success was not promoted to control success; no automatic retry, model fit or native scoring followed.

The [read-only diagnostic](../data/processed/secondary-feature-mapping-control-diagnostic.json), produced by the [audit helper](../scripts/audit-secondary-mapping-control.py), establishes:

- Training protein `htg0000p0000695`, its copied query and the copied query's selected target have the same 385-aa sequence and SHA256 `1fdaf920d0bb4b436db9fa8653ee516cd0f340fb77652b84db5231c19d63c12e`.
- The two focal searches return the same 47 targets in the same order. Forty-six complete target-hit tuples agree.
- The single changed tuple is the original protein's self target: reported identity changes from `1.00` to `0.923`. Score 745, E-value `1.408E-243` and full-length coordinates remain equal.
- The frozen parser's float32 identity/E-value ordering then selects a different family in the external query. Raw-first selection among equal best keys propagates the mismatch into the vector. This is not a dropped candidate or corrupted query sequence.

## Pinned implementation explanation

The installed binary reports commit `8cc5ce367b5638c4306c2d7cfc652dd099a4643f`. The official [alignment implementation](https://github.com/soedinglab/MMseqs2/blob/8cc5ce367b5638c4306c2d7cfc652dd099a4643f/src/alignment/Alignment.cpp#L375) identifies a same-database self hit and explicitly sets its identity to 1. The [matcher](https://github.com/soedinglab/MMseqs2/blob/8cc5ce367b5638c4306c2d7cfc652dd099a4643f/src/alignment/Matcher.cpp#L107) estimates identity from the alignment score in score/coverage mode, while the sequence-identity mode uses aligned identical-residue counts. Both cached source files have verified Git blob hashes, retained in the failure receipt.

The effective v1 search logs use alignment mode 2. The [release-specific official guide](https://github.com/soedinglab/MMseqs2/wiki/Home/67ba9c6637b4b5121a73e5de034dd0c3414d2b81#how-does-mmseqs2-compute-the-sequence-identity) distinguishes the default score-based estimate from mode 3's aligned-residue identity. These documented semantics explain the recorded context difference; the observation is not evidence that MMseqs returned different protein sequences.

## Focal controlled diagnostic

The [probe policy](../software/tipcraft/benchmarks/receptor-resource/mmseqs-identity-mode-probe-policy.json) was hash-frozen before its commands. It selects the same focal record from each existing cached prefilter, preserves the full original query/target databases and their database keys, and compares modes 2 and 3 in training-self and external-copy contexts. All other alignment settings remain fixed. Ten small diagnostic commands complete successfully; no phenotype or native labels are used.

The [probe receipt](../data/processed/mmseqs-identity-mode-probe.json) records:

| Alignment mode | Original target identity, self context | Original target identity, external context | Selected families equal |
| --- | ---: | ---: | --- |
| 2, score/coverage estimate | 1.00 | 0.923 | No |
| 3, aligned-residue identity | 1.00 | 1.00 | Yes |

Mode 3 changes the focal training assignment to match the external assignment; it does not rewrite the failed v1 vocabulary or receipt. This is a causal technical diagnostic on one sequence and its cached candidate sets, not a full copied-genome replay, a new prediction result or permission to count the failed feature task as covered.

## Next justified step and decision

**HOLD biological interaction-model fitting under v1.** Before resuming a production experiment, prospectively specify a consistent training/query identity procedure and validate the complete copied-proteome control, including independent technical copies where practical. A change must receive new versioned policies, fresh destinations and an explicit deviation record; the v1 failure, plan and raw searches remain intact. Source partitions, infection labels, model settings and screening thresholds must not be selected using native outcomes.

No new source or Ando infection score exists. New-method Ando coverage remains **0/10**, WIsH's measured failure remains, and four exact exchange-construct cases continue to abstain. The present failure concerns representation consistency before prediction; it neither establishes a biological failure of all protein-family methods nor validates binding or tail-fibre assembly.
