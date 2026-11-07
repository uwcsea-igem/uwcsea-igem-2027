# Receptor-model inputs and independent evaluation readiness

Updated 2026-11-07. **The accessible receptor resource supports input preparation, while label version and relatedness remain explicit admission gates.** This route uses deposited phage sequences without the declined Picard host-genome collection. It complements the original strain-ranking question; no new prediction coverage is claimed.

## Literature-backed method and endpoint

The [v1 study Methods](https://www.biorxiv.org/content/10.64898/2026.04.02.716166v1.full) describes 13 binary dependency classifiers built from amino-acid 5-mers across phage proteomes. It reports repeated random 10% holdouts for feature generation, selection and classification, and a separate 206-phage development/49-phage prospective cohort. Those statements do not establish sequence-group independence. The final models incorporate the prospective cohort, so they cannot supply an unseen-cohort result. Exact run seeds, commands and a matching inference bundle remain unverified.

This method has a mechanism-relevant endpoint covering both protein and LPS dependencies. Its labels nevertheless derive from genetic-dependence experiments, rather than direct physical binding measurements. A receptor-class prediction alone supplies no validated mapping to the three target strains and no evidence of tail-fibre assembly compatibility.

## Source-label version audit

The official [v1 supplementary tables](https://doi.org/10.6084/m9.figshare.31930314.v1), file 63415071, were downloaded and authenticated: 72,821 bytes, MD5 `cf0b851dc65faf55cff6aa99e444560f`, SHA256 `10de34fd37a7ba1d355b56cdb8b214b690fc71f8991108669dee9d223ab298b7`.

The [frozen comparison policy](../software/tipcraft/benchmarks/receptor-resource/label-source-policy.json) selects all rows and literal fields before comparison. The [executed label-source audit](../data/processed/receptor-label-source-audit.json) finds:

| Comparison of v1 Table S1 and pinned updated browser metadata | Result |
| --- | ---: |
| Source rows / browser rows / unique paired primary IDs | 260 / 260 / 260 |
| Unresolved or ambiguous primary IDs | 0 |
| Differing model-role or susceptibility cells | 0 |
| Differing BW25113 receptor cells | 2 |
| Differing BW25113 LPS-sugar cells | 5 |
| Differing BL21 receptor cells | 5 |
| Distinct primary IDs with a differing annotation | 11 |

These are literal version differences, not adjudicated corrections. The audit preserves both versions and every difference. A future fit must choose and bind one source version before scoring; it cannot silently combine favorable annotations. `Unknown`, `Not assayed`, `Resistant`, `None`, generic `LPS`, and compound receptor annotations remain distinct literal states. An unlisted receptor is not automatically a verified negative class. The supervised target encoding must establish which source evidence supports positive, negative and unknown values for each task.

The [bounded asset-path audit](../data/processed/receptor-model-asset-path-audit.json) inspected the complete 1,349-entry pinned browser repository tree. No non-AlphaFold path matched the declared classifier/feature filename search. This is a path inventory, not proof that training assets are absent from every deposit. Generic GenoPHI code alone does not resolve this study’s exact task-label encoding.

## Authenticated cohort inputs

The [executed cohort inventory](../data/processed/receptor-cohort-container-audit.json) authenticates all 260 cohort GenBank files plus the separately marked catalogue T3 candidate against pinned sizes and Git blob hashes. It records **268 nucleotide records and 48,083 deposited protein translations**, with zero payload errors, parser warnings or CDS missing translation. This is annotation availability, not yet translation agreement or stock identity.

Three containers have multiple records: K20 (three), Mnm1 (five) and TP9 (two). Every record remains separate in the raw FASTA and is mapped to its source container; they are not concatenated or trimmed. Four containers carry ambiguous sequence: JK36, JK42, T4 and TP15. Across them the inventory counts 205 ambiguous nucleotide characters and 71 `X` amino-acid characters in five proteins. The strict k-mer kernel currently rejects such nonstandard protein inputs. Any future ambiguity policy requires an explicit input rule and accounting; these records are not silently cleaned or removed.

The full per-protein metadata remains in an immutable ignored manifest, hash-bound by the compact versioned report. The original full execution report is retained in raw storage; the compact report records this provenance-only refinement rather than claiming the smaller report was the original execution artifact.

## Direct global k-mer projection

The new [k-mer kernel](../software/tipcraft/src/tipcraft/kmer_features.py) provides a label-free implementation candidate with explicit boundaries:

- Generate contiguous 5-mers within each deposited protein; never join proteins.
- Fit the vocabulary using training proteins only. The proposed rare filter requires occurrence in at least two distinct training proteins, which can belong to one genome. Repetition within one protein counts once.
- Merge k-mers with identical training-sample occurrence patterns. Preserve every member of each merged feature and project by their union in new samples.
- Keep vocabulary and feature order fixed during projection. Held-out-only k-mers and held-out occurrences cannot rescue training singletons. Preserve zero vectors and reject missing or nonstandard protein inputs rather than silently cleaning them.

This is a declared implementation adaptation. Its direct global projection has not been certified equivalent to the authors' novel-genome protein-family assignment workflow. It has not been applied to the real cohort or used for feature selection, fitting or scoring. Tests check protein boundaries, protein-level occurrence counting, held-out vocabulary exclusion, merged-member divergence, deterministic ordering and invalid-input rejection.

## Required evidence before scoring

1. Authenticate all cohort containers and inventory deposited translations. Keep failed or incomplete records visible; do not generate missing proteins or conflate catalogue names with experimental stocks.
2. Screen the original target gp17 sequences against cohort DNA with fixed thresholds and query-completion controls. Direct near-match exclusion is only one component: whole-genome/related-group analysis must still establish training and evaluation groups.
3. Resolve source-label version and negative/unknown semantics. Freeze task membership, preprocessing, feature selection, model parameters, grouped splits and matched baselines before a result run. Related parent/mutant phages stay together.
4. Report the receptor endpoint separately. The ten Ando observations, WIsH negative result, four construct abstentions and original strain-ranking decision criterion remain unchanged.

Reproduce the literal label comparison with fresh destinations:

```bash
uv run --project software/tipcraft python scripts/audit-receptor-label-sources.py --workdir data/raw/receptor-resource/label-sources-reproduction-1 --output /tmp/receptor-label-source-audit.json
```

The current disposition is **GO for authenticated input and method preparation; HOLD for fitting or claiming strain-specific predictive performance**.
