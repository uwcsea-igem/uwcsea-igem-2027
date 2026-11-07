# December decision-gate report

Computational evidence update · run date 2026-11-07 · team review pending.

**Decision supported today: NO-GO for using this WIsH model as TipCraft's strain-specific ranking engine; GO for reproducible evidence/evaluation work.** The pipeline now runs an independent, literature-backed sequence scorer and exposes a negative result. Neither the universal socket nor assembly prediction has been validated. This result does not invalidate all other modelling methods or determine wet-lab feasibility.

The [alternative route](22-alternative-route.md) now proceeds without the unavailable Picard host-genome collection. Authenticated receptor-study tables, complete workbook/browser correspondence and two existing target GenBank containers support further input and method validation. This is source-feasibility progress: new independent-method coverage of the ten Ando observations remains **0/10**. The receptor-annotation route requires an additional justified host-context model before it can address strain rankings; its accuracy would not establish assembly compatibility.

The [host-encoder capacity test](26-host-encoder-capacity.md) additionally supports **NO-GO for the proposed two-reference, merged5-mer/union host representation as the sole host input to a paired ranking engine**. Its three original target-host vectors are identical in both the deposited-annotation pilot and a separately fixed uniform-Prodigal control, so a downstream model cannot distinguish those hosts through those inputs. This is a label-free representation result, not a new infection score, and does not reject every encoder or all GenoPHI methods.

A separately frozen [family-conditioned binary host pilot](29-family-host-capacity.md) passes its source self-projection controls and distinguishes BL21 from the K-12 references, but BW25113 and MG1655 still have identical vectors. That particular representation therefore cannot distinguish those two hosts as the sole host input. No classifier was fitted and no threshold search followed the result; independent-method coverage remains 0/10. This is a specific representation boundary, not a universal failure of family-aware methods.

The [endpoint correction](30-collision-endpoint-correction.md) clarifies that this remaining family-vector collision does **not** rule out the native ranking task: BW25113 and MG1655 share labels for both native phages. All-three vector distinction was an overstrict admission requirement. T3 is the sole informative native ranking query, so a host-only comparator remains necessary before attributing any future ranking signal to phage-specific input. Source/training admission remains HOLD; the measured global-union failure is unchanged.

The [independent source-host audit](31-source-host-contrast-readiness.md) now identifies an additional validation limitation: only six of 49 candidate testing phages have complete literal two-host states, and both informative ordering comparisons favor BW25113. A fixed host preference can satisfy all those comparisons. This subset alone therefore cannot validate a phage-dependent ranking advantage. Missing states remain unmeasured; no new predictor or accuracy result is claimed.

The subsequent [BASEL–ECOR source audit](32-secondary-contrast-source-selection.md) identifies 167 host pairs with opposite observed directions across phages in a different, published 52-phage/25-host panel. This supports continued admission work for a secondary validation source. Genomic exclusions, splits, stock correspondence and spot-clearance endpoint limits remain unresolved; no new model or original Ando coverage follows.

All 25 secondary-panel [ECOR WGS inputs are now authenticated](33-ecor-wgs-input-authentication.md), retaining 3,477 versioned contigs with a complete offline replay of 100 raw artifacts. This removes a reference-input prerequisite without using the skipped Picard collection. It does not establish genetic independence, stock equivalence or prediction; grouping and phage accession checks remain pending.

## Method and fixed inputs

We used the [authors' WIsH executable](https://github.com/soedinglab/WIsH/tree/10e3f11779e926fb6dbdb20197c248367a75887f), version 1.1, with its default order 8 and alpha 16, built in Release mode with Apple clang 21, single-threaded without OpenMP. [Galiez et al., 2017, Methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC5870724/) motivates fitting one Markov sequence model per candidate host and selecting the host with maximum likelihood. No known interaction labels or donor mapping are used to fit or score. The original paper evaluates taxonomic host prediction, including genus-level performance. This run is an exploratory strain-level stress test beyond that demonstrated task. It does not score gp17 structure or assembly.

The method, candidate set, parameters and endpoints were recorded before the first score run in [wish-inputs.json](../software/tipcraft/benchmarks/retro-t7-t3/wish-inputs.json). Labels were already visible when the method was selected: this is exploratory retrospective evaluation, not a preregistered or blinded test. Parameters were not optimized after seeing these scores.

| Input | Public accession | Length (bp) |
| --- | --- | ---: |
| T7 | NC_001604.1 | 39,937 |
| T3 | NC_003298.1 | 38,208 |
| BL21, TaKaRa reference | CP010816.1 | 4,528,118 |
| BW25113 | CP009273.1 | 4,631,469 |
| MG1655 | U00096.3 | 4,641,652 |

The manifest records source URLs, FASTA headers, sequence/file hashes and retrieval date. The public references are not demonstrated to be identical to the experimental stocks. In particular, BL21 was not silently replaced with BL21(DE3). Four exchange cases remain abstentions: no verified complete-genome accession was mapped to those exact constructs, and we did not invent or construct a hybrid genome to fill that gap. The [follow-up source audit](08-construct-source-audit.md) records the bounded search and inaccessible supplements; it does not establish that no public record exists.

## Actual independent scores

Higher mean conditional log likelihood ranks higher. These are the unmodified upstream output values, rounded by WIsH to six significant digits. They are not binding energies, calibrated infection probabilities or assembly scores. No Gaussian null distribution was fitted, so no WIsH p-values are available. Compare host scores within each phage; do not pool their magnitudes across phages.

| Phage | Candidate host | Mean log likelihood | Host rank | Published plaque label |
| --- | --- | ---: | ---: | --- |
| T7 | BW25113 | -1.41957 | 1 | Detected |
| T7 | MG1655 | -1.41971 | 2 | Detected |
| T7 | BL21 | -1.42012 | 3 | Detected |
| T3 | BW25113 | -1.41117 | 1 | Not detected |
| T3 | MG1655 | -1.41128 | 2 | Not detected |
| T3 | BL21 | -1.41199 | 3 | Detected |

The qualitative labels are curated from the native-phage observations in [Ando et al., 2015, host-range section and Fig. 3](https://pmc.ncbi.nlm.nih.gov/articles/PMC4785837/). Below-detection labels are conditional on that experiment, not absolute claims of noninfectivity.

Primary endpoint: a top-ranked host is a documented positive for **1 of 2 phage queries (50%)**. T7's three candidates are all positive, so that query cannot distinguish useful from useless rankings; a uniformly random top host would also succeed for T7. Across this fixed set, the random-top-host expected rate is 2/3, not 1/2.

Secondary endpoint: T3 places its one documented positive below both negatives, giving **0/2 positive–negative concordant comparisons (0%)**. T7 has no negative among its candidates, so its concordance is undefined. The always-tied score reference has expected concordance 50%.

For compatibility with the previous binary baseline only, a forced-choice sensitivity analysis treats the top host as positive and all others as negative. On the same **six covered observations**, it gives accuracy **33.3%** and balanced accuracy **37.5%**, versus **66.7%** and **50%** for always-positive. WIsH does not natively predict a complete susceptibility matrix: these binary numbers are an explicitly artificial adaptation, not the primary validation. The earlier 60% baseline used all ten observations and must not be compared to this six-case result.

## Leakage and integrity checks

- Scoring accepts only the genome manifest, five FASTA files and an external executable. It has no phenotype CSV input. Headers are stripped to a neutral value before either run. Evaluation reads labels only after the scores have been saved.
- We actually reran upstream model construction and prediction after replacing every host and phage filename with anonymous IDs, then mapped them back. All six scores were identical at upstream output precision. The evaluator checks named/anonymous run flags and matching method, executable, manifest and parameter provenance, and rejects reuse of a named run as the anonymous control. Names and metadata did not supply the result. With `--manifest`, the report also binds both score JSON files to the actual frozen input manifest and their adjacent upstream likelihood matrices. It rejects inconsistent exclusions, candidate sets, parameters and score values. A run declaring label access cannot receive a passing audit. Omitting the optional manifest preserves the older command interface but yields only metadata-level checks, never verified artifacts.
- We froze and hashed the scores, then exhaustively moved T3's single positive label among its three candidate hosts. Evaluation concordance became 0, 1 and 0.5 while the score hash remained unchanged. The upper-tail fraction for the observed result is 3/3; this tiny permutation check supplies no evidence of above-random performance. It is not a powered significance test.
- The five reference sequences have distinct sequence hashes. Candidate-host genomes are intentionally used as unsupervised model inputs; they are not phenotype training examples. No positive/negative association database, pretrained interaction-model weights, train/test row split, donor phenotype lookup or inferred chimera was used.
- Residual bias remains: method and cases were selected with prior knowledge of outcomes; two native phages and closely related hosts are tiny correlated samples; public-stock identity is uncertain. The report and qualitative label mapping received an independent AI-agent review against the primary sources; independent human figure review remains pending. No generic claim of “leakage-free validation” is justified.

## December gate disposition and next evidence

| Question | Current evidence | Disposition |
| --- | --- | --- |
| Does a real label-independent scorer run? | Six nonconstant sequence scores, raw matrix, provenance and repeatable commands | Yes, computational deliverable exists |
| Does it distinguish the decisive native T3 cases? | T3 positive ranks last | NO-GO for this method as the strain-ranking engine |
| Does it recover the published engineered exchanges? | Four explicit abstentions | Unanswered; current 60% coverage fails the proposed complete-coverage gate |
| Does it predict gp17 assembly or arbitrary binder grafting? | No direct assay or calibrated structural model | HOLD; no assembly success rate claimed |
| Can wet-lab work start? | Exact materials, access, supervisor and budget remain unverified | Separate readiness decision, unchanged by these scores |

Retain the failed method as a transparent comparator. Any next method needs a documented relationship to receptor-level specificity, independently sourced calibration examples, exact construct provenance and a frozen held-out evaluation. Do not pick a donor-identity lookup or retune this method until ten familiar labels look right. Gate dates remain proposals; this report is evidence for the team, not an autonomous team decision.

## Artifacts and reproduction

[Benchmark commands](../software/tipcraft/benchmarks/retro-t7-t3/README.md) · [machine-readable report](../data/processed/retro-t7-t3-wish.json) · [raw named-run likelihood matrix](../data/processed/wish/named/llikelihood.matrix) · [anonymous-run scores](../data/processed/wish/anonymous/scores.json) · [build provenance](../data/processed/wish/build-provenance.json).

Raw reference genomes, source archive, compiled GPL tool and model files remain under ignored `data/raw/wish/`. Only the original MIT wrapper, small numeric outputs and provenance are versioned. WIsH retains GPL-3.0; its source or executable is not redistributed as MIT. No wet-lab procedures, new designs or orderable DNA were produced by this evaluation.

## Next-method review, 2026-11-07

The [follow-up feasibility review](09-next-method-feasibility.md) identifies GenoPHI as a relevant strain-level candidate, rather than transferring a genus-level predictor or treating rough LPS as a protein receptor. A pinned [public-input inventory](../data/processed/genophi-input-audit.json) confirms external E. coli data exist, including a host named BL21, but does not establish target/homology exclusion or an authenticated inference bundle. **HOLD for a new independent run**, with zero new-method coverage and no new performance estimate. This is an input-readiness finding, not a GenoPHI failure. WIsH's measured six scores, NO-GO and four exchange abstentions remain unchanged. The [next-run contract](../software/tipcraft/benchmarks/retro-t7-t3/next-method-contract.json) records evaluation rules and explicitly unresolved prerequisites; it is not a completed executable preregistration.

The subsequent [relatedness audit](15-relatedness-audit.md) actually computes the public phage comparisons, identifies a hidden near-case (`427_P2`) for provisional exclusion, and reconstructs related-host groups. Direct execution of the original R function verifies the directional calculation on all 4,851 frozen sequence pairs. The [annotation follow-up](16-annotation-feasibility.md) verifies DNA correspondence for all 2,980 deposited example proteins and fixes uniform reannotation as a declared adaptation. Full host ANI exclusions and the complete feature/split/training specification remain open; these input-curation results do not change the host-prediction or assembly gate dispositions.

The [feature-representation pilot](17-feature-representation-feasibility.md) subsequently learns a vocabulary from 93 source phages without target proteins, yielding different nonzero T7/T3 vectors and an exact copied-proteome mapping control. Neither reference gp17 is assigned a feature. This supports the feasibility of the native whole-genome feature path while withholding any claim of recognition-head scoring. No interaction model or new host score has been produced; coverage and gate dispositions remain unchanged.

The [fitting-component checks](18-training-protocol-and-component.md) now execute RFE/CatBoost on synthetic data and verify training-only interfaces, host-only controls and frozen-bundle authentication. They do not supply biological performance or complete the real-data training specification. The [compatibility catalogue](19-compatibility-evidence-catalogue.md) adds eight published binding/infection observations, with all direct assembly labels unknown; it cannot yet calibrate an assembly scorer.

The [source-admission stage](20-source-admission-and-partitions.md) now authenticates 8/403 host files and returns HOLD for 395 missing files and the unrun host ANI audit; it writes no final roster or partitions. The [receptor-source check](21-receptor-dependence-and-source-check.md) distinguishes the team's intended native-LPS and protein-targeting questions and corrects categorical knockout-control assumptions using published T7 data. Neither follow-up adds target predictions or changes the measured WIsH result.

The [grouped source-phage feature run](17-feature-representation-feasibility.md#executed-grouped-source-phage-preprocessing) completes three separately learned vocabularies with passing copied-proteome controls. Nine of 33 held-out phages in one fold have no assigned proteins and identical zero vectors; the rows and original settings are retained. This identifies a representation limitation to carry into attribution/evaluation. Host features and interaction models remain unrun, so the independent-prediction gate remains HOLD with 0/10 new-method coverage.
