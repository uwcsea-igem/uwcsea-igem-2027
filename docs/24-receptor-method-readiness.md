# Receptor-model inputs and independent evaluation readiness

Updated 2026-11-07. **The accessible receptor resource supports input preparation, while label version and relatedness remain explicit admission gates.** This route uses deposited phage sequences without the declined Picard host-genome collection. It complements the original strain-ranking question; no new prediction coverage is claimed.

## Literature-backed method and endpoint

The [v1 study Methods](https://www.biorxiv.org/content/10.64898/2026.04.02.716166v1.full) describes 13 binary dependency classifiers built from amino-acid 5-mers across phage proteomes. It reports repeated random 10% holdouts for feature generation, selection and classification, and a separate 206-phage development/49-phage prospective cohort. Those statements do not establish sequence-group independence. The final models incorporate the prospective cohort, so they cannot supply an unseen-cohort result. Exact run seeds, commands and a matching inference bundle remain unverified.

This method has a mechanism-relevant endpoint covering both protein and LPS dependencies. Its labels nevertheless derive from genetic-dependence experiments, rather than direct physical binding measurements. A receptor-class prediction alone supplies no validated mapping to the three target strains and no evidence of tail-fibre assembly compatibility.

Target outcomes and source annotations were already visible during method development. Freezing a future run prevents result-driven changes to that run; it does not turn this work into a blinded study.

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

Three containers have multiple records: K20 (three), Mnm1 (five) and TP9 (two). Every record remains separate in the raw FASTA and is mapped to its source container; they are not concatenated or trimmed. Four containers carry ambiguous sequence: JK36, JK42, T4 and TP15. Across them the inventory counts 205 ambiguous nucleotide characters and 71 `X` amino-acid characters in five proteins. The default strict k-mer profile rejects nonstandard protein inputs. The separately frozen deposited-sequence profile below now gives an explicit window-level rule, without cleaning or removing those source samples.

The full per-protein metadata remains in an immutable ignored manifest, hash-bound by the compact versioned report. The original full execution report is retained in raw storage; the compact report records this provenance-only refinement rather than claiming the smaller report was the original execution artifact.

## Direct global k-mer projection

The new [k-mer kernel](../software/tipcraft/src/tipcraft/kmer_features.py) provides a label-free implementation candidate with explicit boundaries:

- Generate contiguous 5-mers within each deposited protein; never join proteins.
- Fit the vocabulary using training proteins only. The proposed rare filter requires occurrence in at least two distinct training proteins, which can belong to one genome. Repetition within one protein counts once.
- Merge k-mers with identical training-sample occurrence patterns. Preserve every member of each merged feature and project by their union in new samples.
- Keep vocabulary and feature order fixed during projection. Held-out-only k-mers and held-out occurrences cannot rescue training singletons. Preserve zero vectors; reject missing proteins and characters outside the chosen profile rather than silently cleaning them.

The [candidate representation policy](../software/tipcraft/benchmarks/receptor-resource/direct-kmer-representation-policy.json) now fixes an explicit `deposited21_skip_X` profile after source-character inspection and before any biological feature fit. Both host protein files contain three literal U residues. The profile retains U as distinct from unknown X, skips only 5-mer windows containing X, and counts repeated k-mers once per original protein even across X-separated regions. It never deletes X, joins its neighbors, creates independent fragment proteins or replaces deposited first residues. Projection reports skipped-window and X/U counts for every sample. The original `standard20` default remains strict. Synthetic tests cover these differences; no real-cohort feature vocabulary has yet been fitted.

This is a declared implementation adaptation. Its direct global projection has not been certified equivalent to the authors' novel-genome protein-family assignment workflow. It has not been applied to the real cohort or used for feature selection, fitting or scoring. Tests check protein boundaries, protein-level occurrence counting, held-out vocabulary exclusion, merged-member divergence, deterministic ordering and invalid-input rejection.

## Executed annotation and near-homolog checks

The [strict translation-consistency report](../data/processed/receptor-cds-translation-check.json) covers all 48,083 deposited translations under a policy fixed before checking:

| Status under the declared full-CDS translation convention | Entries |
| --- | ---: |
| Exact deposited/computed string match | 44,212 |
| String mismatch | 3,666 |
| Full-CDS validation failed | 199 |
| Unverified because coding DNA is ambiguous | 6 |

A separate [post hoc descriptive follow-up](../data/processed/receptor-cds-difference-description.json) finds that **all 3,666 mismatches have the same length and differ only at the first residue**. The complete-CDS translation convention gives initial M; the corresponding deposited residues are V (1,994), L (1,669) or I (three). This pattern is consistent with alternative-initiation representation conventions, rather than broad internal protein/DNA disagreement; it does not independently establish how the annotations were generated. The original mismatch classifications and sequences remain unchanged. Of the 199 full-CDS validation failures, 170 lack the expected terminal stop codon and 29 lack a recognized start codon under this check. They are not infection or receptor failures. A future representation policy must address these annotation conventions and unresolved entries explicitly.

The [executed DNA-based gp17 screen](../data/processed/receptor-gp17-overlap.json) preserves all 268 cohort/catalogue nucleotide records and appends two separate original target-reference controls. The database count of 270, both query-completion records and both reference self-controls passed. With the frozen single-alignment rule (at least 90% identity and 80% full-query coverage), **Bas66, T7 and the catalogue T3 candidate** meet direct near-hit exclusion; each has a 100%-identity, full-query-coverage hit. Bas66 and T7 are source training members. Bas64 additionally meets the 40%/80% broad diagnostic rule, which is not an automatic near-hit exclusion.

The 257 candidates without a reported hit meeting either rule are **not** cleared of whole-genome relatedness, distant homology or group overlap. No final training roster is admitted by this direct screen. The annotation-consistency findings do not change the DNA screen or cause samples to disappear from it.

## Host-context coverage cannot be inferred from propagation metadata

The [executed host-context audit](../data/processed/receptor-host-context-audit.json) retains all 260 metadata rows and finds explicit susceptibility fields for **BW25113 and BL21 only**, with no corresponding MG1655 field. Production metadata contains `MG1655_RM` (153 rows) and `MG1655_RM wbbL(+)` (33 rows); these are not additional susceptibility measurements and are not used as predictor features or positive labels.

Official v1 Table S2 lists `MG1655ΔRM` with a Maffei et al. 2021 reference. That [primary study, Bacterial strains and strain construction](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3001424) describes a restriction-system-deficient MG1655 derivative and separately altered surface-glycan backgrounds. This establishes that derivative context matters; it does not establish stock equivalence between differently formatted names in the two resources. Neither derivative name is silently equated with the original Ando MG1655 case. Differences in intracellular defence also reinforce why receptor dependence and functional infection remain separate endpoints.

The BL21 susceptibility column includes 96 `Not assayed` values and one literal `Resistant`, alongside `Yes` and `No`. The audit preserves all states without binary conversion. A future known-host diagnostic would require its own justified label semantics and would still not supply a validated score for the original third host. The existing six-native-observation criterion and ten-observation coverage denominator therefore remain unchanged.

The [host-asset candidate audit](../data/processed/receptor-host-asset-candidates.json) identifies six pinned nucleotide/protein/annotation assets in the author browser, independently of Picard. Two small browser configurations are byte-authenticated: BW25113 loads `Keio`, while BL21 names `CP001509`. The [original genome publication](https://pubmed.ncbi.nlm.nih.gov/19786035/) identifies CP001509 as BL21(DE3). At the candidate-audit stage the payloads were uninspected. The subsequent authentication below establishes file identity and declared sequence headers, while experimental-stock correspondence remains unverified. The original BL21 TaKaRa reference is retained.

This leaves a concrete paired host-feature route to investigate after payload authentication and source correspondence. Missing MG1655 training annotations do not mathematically prevent an explicitly specified unseen-host model from producing a score; they do prevent treating production metadata as its ground truth. Training on two source host backgrounds would require clear limits on unseen-host validation and comparison with host-only/prevalence baselines. No such paired model has been specified or fitted here.

The [executed host-payload audit](../data/processed/receptor-host-payload-audit.json) subsequently authenticated all six listed files against pinned Git blob hashes and inventoried them without sequence editing:

| Browser source directory | Existing nucleotide record | Length | Deposited protein records | Gene-table rows |
| --- | --- | ---: | ---: | ---: |
| BL21 | CP001509 | 4,558,953 | 4,195 | 4,303 |
| Keio | 7023 | 4,639,675 | 4,146 | 4,610 |

Both are single-record nucleotide files, and both sequence hashes differ from the corresponding original public Ando-case host references. The protein and gene-table counts describe different source objects and are not automatically a one-to-one CDS match. Each source protein file contains three U residues, retained explicitly by the candidate profile above. These author-browser reference representations are now available for a source-correspondence check, not automatically admitted experimental genomes. No source susceptibility label has been reassigned to another strain, and no original target reference has been replaced.

The preprint's cited library source, [Mutalik et al. 2020, Methods](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3000877), explicitly distinguishes its BL21 RB-TnSeq library from a BL21 Dub-seq library mapped to a BL21-DE3 reference and assayed in a BL21DE3C43 background. This supports treating reference mapping and assay stock as separate provenance fields. It does not establish that the 2026 susceptibility annotations or RB-TnSeq experiments used the Dub-seq background. The modalities and their source correspondences must not be merged by their shared BL21 name.

```bash
uv run --project software/tipcraft python scripts/audit-receptor-host-payloads.py --cache data/raw/receptor-resource/host-payloads-v1 --output /tmp/receptor-host-payload-audit.json
```


## Required evidence before scoring

1. Container authentication and translation inventory/checking are complete as reported above. Freeze how unresolved annotations, ambiguous residues and multiple records will be represented before feature construction; do not silently correct or remove them.
2. The direct target-gp17 screen is complete. Whole-genome/related-group analysis must still establish training and evaluation groups, including parent/mutant and nonmodel bridge records.
3. Resolve source-label version, negative/unknown semantics and host-context applicability. Freeze task membership, preprocessing, feature selection, model parameters, grouped splits and matched baselines before a result run. Related parent/mutant phages stay together.
4. Report the receptor endpoint separately. The ten Ando observations, WIsH negative result, four construct abstentions and original strain-ranking decision criterion remain unchanged.

Reproduce the literal label comparison with fresh destinations:

```bash
uv run --project software/tipcraft python scripts/audit-receptor-label-sources.py --workdir data/raw/receptor-resource/label-sources-reproduction-1 --output /tmp/receptor-label-source-audit.json
```

The current disposition is **GO for authenticated input and method preparation; HOLD for fitting or claiming strain-specific predictive performance**.

Reproduce the two input checks with fresh output locations:

```bash
uv run --project software/tipcraft --with biopython==1.86 python scripts/check-receptor-cds-translations.py --workdir data/raw/receptor-resource/cds-reproduction-1 --output /tmp/receptor-cds-translation-check.json
uv run --project software/tipcraft python scripts/audit-receptor-gp17-overlap.py --workdir data/raw/receptor-resource/gp17-reproduction-1 --output /tmp/receptor-gp17-overlap.json
```

The dataset-bound descriptive replay uses the byte-identical [helper](../software/tipcraft/benchmarks/receptor-resource/cds-diagnostic-v1/describe.py) and adjacent [post hoc plan](../software/tipcraft/benchmarks/receptor-resource/cds-diagnostic-v1/frozen-plan.json). Copy both to a fresh raw working directory and run the helper there from the repository root, with Biopython 1.86 and the authenticated baseline caches available. It creates `diagnostic.json` beside itself and refuses an existing output; it never updates the baseline check.
