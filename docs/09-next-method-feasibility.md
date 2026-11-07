# Next-method feasibility review

Reviewed 2026-11-07; bounded primary-source and input audit, before any new model scores. **HOLD for another independent T7/T3 evaluation; prioritize GenoPHI for input preparation.** This is an input-readiness decision, not a measured failure of GenoPHI, a claim that no suitable method exists, or team approval. The previous WIsH NO-GO remains in force.

## Biological and methodological fit

[González-García et al. (2015), Abstract and Results](https://pmc.ncbi.nlm.nih.gov/articles/PMC4400320/) demonstrate T7 fiber recognition of rough LPS and LPS-triggered DNA ejection in vitro. A protein–protein docking score therefore has no automatic physical interpretation for native T7. Biosynthesis-enzyme sequences can be statistical proxies for surface glycans; those enzymes are not the extracellular receptor. Receptor access, adsorption, post-entry infection and fiber assembly remain separate endpoints. None of the reviewed host predictors measures assembly.

| Candidate | Demonstrated task / representation | Decision for this benchmark |
| --- | --- | --- |
| **GenoPHI**, Noonan et al. 2026 | Strain-level genome-pair prediction including E. coli; protein-family features, CatBoost and evaluations involving unseen phages and hosts. [Primary paper](https://doi.org/10.1038/s41564-026-02482-5), Methods: model configurations and novel-genome assignment | **Preferred conditional candidate.** Relevant empirical resolution; no requirement to misrepresent LPS as a protein target. This supports testing, not a T7/T3-specific performance claim. Input and overlap prerequisites below are unmet. |
| Gaborieau et al. 2024 per-phage classifiers | E. coli host genomic features, including adsorption factors. Each classifier is specific to a characterized phage and has no phage sequence features. [Paper](https://doi.org/10.1038/s41564-024-01832-5); [owning implementation](https://github.com/mdmparis/coli_phage_interactions_2023) | Cannot instantiate a novel T3 classifier just from its genome. Original [Guelin metadata](https://github.com/mdmparis/coli_phage_interactions_2023/blob/main/data/genomics/phages/guelin_collection.csv) contains T7_Portugal (39,778 bp), whose identity is not established as our 39,937-bp T7 reference or Ando stock. Transferring a related phage's phenotype is not independent prediction. |
| PhageHostLearn, Boeckaerts et al. 2024 | Klebsiella RBP/K-locus capsule features, grouped strain evaluation. Spot-test labels do not necessarily establish productive infection. [Primary paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC11111740/), Methods: interaction data and LOGOCV | Useful evaluation precedent; no demonstrated transfer of its trained capsule model to rough E. coli LPS. Do not import its published accuracy. |
| PhageMind, Shen et al. 2026 | Tail fibers, O-antigen biosynthesis clusters, known infection graph and labelled support examples. [Primary paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC13340165/), §§2.1–2.3; [code](https://github.com/YangSH-ac/PhageMind) | Conditional lead. Rough-LPS applicability and target-excluded training/graph construction remain unverified. Edge masking and random pair partitions alone do not establish independent construct/publication/homology groups. |
| CHERRY / iPHoP | CHERRY describes species-level training and lack of strain optimization; iPHoP targets genus assignment. [CHERRY Discussion](https://pmc.ncbi.nlm.nih.gov/articles/PMC9487644/); [iPHoP](https://pmc.ncbi.nlm.nih.gov/articles/PMC10155999/) | Do not repeat a taxonomic-resolution stress test as if it were a supported strain-specific method. |

## What was actually audited for GenoPHI

The [authors' source](https://github.com/Noonanav/GenoPHI/tree/f7e30abc1b0bbba1ae76b8f58a3363bbbc5cd54f) is pinned at `f7e30abc1b0bbba1ae76b8f58a3363bbbc5cd54f`. We inspected its recursive file tree, three E. coli interaction tables, inference implementation and the [published processed-data archive](https://zenodo.org/records/21696332). The [machine-readable inventory](../data/processed/genophi-input-audit.json) records URLs, hashes and archive contents.

Actual counts: **37,788 pairs / 402 hosts / 94 phages**, a **16,638-pair / 177-host** subset, and a **625-pair / 25 × 25** demo. Only 25 host and 25 phage demo proteomes are bundled in the inspected repository tree. The full and subset tables contain a host named **BL21**. No T7/T3 name substring matched their phage IDs; that does **not** establish sequence or homology exclusion. Original collection metadata and the filtered modelling tables must not be conflated.

No usual model-file extensions were found in the inspected tree or archive; the live releases endpoint also returned an empty list during review. These bounded observations do not exclude other deposits, nonstandard formats or future releases. The [assignment workflow](https://github.com/Noonanav/GenoPHI/blob/f7e30abc1b0bbba1ae76b8f58a3363bbbc5cd54f/genophi/workflows/assign_predict_workflow.py) requires trained models, reference clusters and a feature map. A phenotype table or SHAP summary is not that inference bundle.

**Consequences:** no authenticated ready-to-use model/reference bundle has been established here. Retraining is a possible route, but needs accession/proteome-to-label mapping and a sequence-based exclusion audit first. BL21's presence is a concrete overlap warning: renaming it cannot make it an unseen host. Demo-only retraining would be a different, much smaller experiment and is not silently substituted for the published method. We have not established that the remaining public data are unavailable.

The audit downloaded public tables containing labels; it is **not label-isolated scoring**. It did not train, choose a cutoff using target outcomes, read the archive's prediction values, or generate host predictions. New-method coverage is **0/10**, with performance **not evaluated**. WIsH remains **6/10**; these are separate method-specific denominators.

### Training-source mapping follow-up

The [source-map artifact](../data/processed/genophi-training-source-map.json) now connects the modelling table to the original collection's pinned repository and versioned host deposit. All **94 phage IDs** have exact-name genome-file metadata matches. Of **402 hosts**, **265** match directly and **137** have unique punctuation-normalized candidates. The latter are recorded for review, not accepted as authenticated aliases. The deposit lists 403 host assemblies (about 2.06 GB); public genome sources therefore exist beyond the demo. Sequence download/verification, annotation correspondence and relatedness exclusions are the next input tasks. These metadata matches do not authenticate experimental-stock identity or supply new predictions.

Reproduce after the first audit populates its matrix cache:

```bash
uv run --project software/tipcraft python scripts/map-genophi-training-sources.py --cache data/raw/genophi-training-map --matrix data/raw/genophi-audit/data/interaction_matrices/ecoli_interaction_matrix.csv --output /tmp/genophi-training-source-map.json
```

Add `--offline` for cache replay. The script authenticates the original collection's Git tree and pinned modelling-table hash, and checks a frozen canonical hash of the reviewed host-file inventory. Deposition metadata includes file IDs and checksums for later sequence retrieval. The local tree cache came from the GitHub connector, with the same API-rate-limit consideration as the first audit.

The [genome-preparation follow-up](13-genome-input-preparation.md) now records authenticated downloads, a portable retrieval kit and the frozen host-triage policy. The [label-provenance audit](14-training-label-provenance.md) cross-checks every shared nonblank outcome and excludes 156 source-empty pairs from future training pending clarification. These are input-curation findings, not new host-prediction results.

The [actual relatedness audit](15-relatedness-audit.md) identifies a near-target source genome, `427_P2`, in the training table despite the absence of a T7/T3 name. It also reconstructs 301 related-host groups. Direct execution of the original R function verifies the directional calculation on all frozen pair data; an unrelated publication-control name/accession conflict remains unresolved. Final exclusions still await the complete host audit and training manifest. These are input-curation results, not host-prediction validation.

The [annotation follow-up](16-annotation-feasibility.md) runs a frozen gene-calling procedure on 110 inputs and recovers both reference gp17 sequences. The regenerated sets differ from all 25 pinned example proteomes, while all 2,980 deposited proteins translate exactly from their declared intervals in the authenticated DNA. A separate component policy now fixes uniform reannotation as an explicit adaptation for every retained training, validation and target genome. This settles the annotation choice without claiming exact upstream replication; full clustering, splits and modelling parameters remain to be frozen before fitting.

The [train-only feature pilot](17-feature-representation-feasibility.md) now builds a 297-feature vocabulary from 93 eligible source phages and assigns different nonzero vectors to T7/T3. Both reference gp17 proteins receive no assignment under those settings. This establishes whole-genome representation capacity without recognition-head coverage; it supplies no host scores. GenoPHI remains a candidate for native-genome retrospective evaluation, not an established scorer of interchangeable heads.

The [training-flow review and fitting component](18-training-protocol-and-component.md) now distinguish reusable fitting primitives from upstream wrapper resplits that do not preserve this project's untouched outer test. The new component actually fits RFE/CatBoost on synthetic data, with no evaluation-set argument and externally anchored bundle scoring. These software checks do not authenticate the biological training membership or complete the production protocol; full source data and fold-specific preparation remain required.

The [source-admission stage](20-source-admission-and-partitions.md) now connects all 402 host IDs through the original source records to exact deposit IDs, derives complete phage groups and defines deterministic crossed host/phage blocks. Its actual run remains HOLD with 395 host files missing and no host ANI audit; no final training roster or partitions were written. This later identifier crosswalk refines the earlier 137 candidate-name status without asserting stock identity.

## Fixed next-run contract and reopening conditions

[The evaluation contract](../software/tipcraft/benchmarks/retro-t7-t3/next-method-contract.json) fixes the candidate, target references, exclusions, leakage rules, comparison subset and interpretation before any new scores. It is explicitly **not an executable preregistration**: training inputs, exact model parameters and inference artifacts are unresolved. Freeze and hash the training specification (inputs, exclusions, splits, versions, parameters and seeds) **before fitting**. Then hash the resulting learned clusters, feature mappings and weights in a trained-bundle manifest **before target inference**. An already trained public bundle instead requires an authenticated training specification and overlap audit before inference. Retrospective declarations cannot substitute for either stage.

To reopen HOLD, supply either an authenticated trained bundle with auditable training membership, or a reproducible retraining dataset with mapped genomes and labels. The audit must exclude target stocks, their related groups and target publication records from supervised training and calibration; fit all preprocessing and feature selection only within training partitions. Record retained/excluded IDs and reasons, not just an assertion of independence. Freeze exact software/dependency versions, parameters, seeds, feature mappings and model hashes. Confirm distinct input feature vectors can represent the intended contrasts; report collisions without changing the representation after viewing target scores.

The primary diagnostic is within-phage positive–negative concordance, compared with tied/random rankings and WIsH **on identical covered cases**. The current informative query is T3; T7's all-positive panel cannot establish discrimination. With just three possible placements of T3's single positive, even perfect ordering has an exact upper-tail fraction of 1/3: it cannot establish reliable generalization. An apparent win is a reason to seek an independently sourced grouped holdout, not to claim a validated engine.

Four exchange observations retain the [existing source abstentions](08-construct-source-audit.md). Missing exact genomes block a whole-genome scorer for those cases; this does not impose a universal full-genome requirement on every possible future method. A fragment-based method would need its own authenticated inputs and pre-score contract. Host prediction cannot reopen the assembly gate.

## Reproduce the input audit

From the repository root, using Python 3.11+ and network access:

```bash
uv run --project software/tipcraft python scripts/audit-genophi-sources.py --cache data/raw/genophi-audit --output /tmp/genophi-input-audit.json
uv run --project software/tipcraft python scripts/audit-genophi-sources.py --cache data/raw/genophi-audit --output /tmp/genophi-input-audit-offline.json --offline
diff -u /tmp/genophi-input-audit.json /tmp/genophi-input-audit-offline.json
```

Raw third-party material stays in ignored `data/raw/genophi-audit/`. CSV blobs are checked against the pinned Git tree; tree hashes are reconstructed and bound to its recorded root; the Zenodo ZIP is checked against the reviewed deposit checksum. The JSON hashes describe cached artifacts, so GitHub response formatting can change the tree-response byte hash without changing its authenticated contents. The local tree response was obtained through the GitHub connector after unauthenticated API rate limiting; a fresh CLI fetch may require waiting for that limit to reset. No credentials are embedded. Offline replay requires the populated cache; this new audit is not a network-free CI check.

This review combines primary-source reading and an independent AI-agent applicability review; human figure/stock-identity review remains pending. No new sequence, design, orderable DNA, binding score or assembly claim was produced.
