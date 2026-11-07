# Independent evaluation without the Picard genome download

Updated 2026-11-07. **The active route is now independent method and evidence selection using accessible public resources. Full Picard-genome retrieval is no longer a prerequisite for project progress.** The repository owner requested this computational rerouting; it does not change the team's biological scope or establish a new prediction result.

## Scientific objective and execution boundary

The central question remains whether an independently specified method supplies useful strain-specific information beyond matched baselines. Preserve the ten Ando observations, WIsH's measured negative result and four construct-source abstentions. The GenoPHI preparation and grouped representation findings remain evidence, while its full-data retraining is optional and on HOLD.

No further Picard whole-host download is scheduled. Alternative archives must not become an unannounced substitution for those inputs. New methods require their own source, overlap and pre-score contracts; the existing GenoPHI contract is retained as a historical, method-specific artifact. Known target outcomes remain excluded from training, calibration and support examples.

## Candidate triage completed

| Resource | Verified opportunity | Next admission requirement |
| --- | --- | --- |
| **PBIP**, Ma et al. 2025 | Official public folders list checkpoints and embeddings, providing an asset-audit route independent of the failed download. | Its strain model was developed on 104 phages and 120 *K. pneumoniae* strains. Transfer to *E. coli* requires evidence; it is not an admitted T7/T3 predictor. Reconstruct and authenticate preprocessing, scalers, sequence membership and checkpoints before scoring. |
| **PhageMind**, Shen et al. 2026 | The linked archive directory lists five checkpoints and processed features. Metadata inspection establishes availability, not model authenticity or target exclusion. | Resolve graph/support-label dependence, training overlap and rough-LPS applicability. Its pinned decoder has the fixed-graph ranking limitation below. |
| **Moriniere et al. 2026 receptor resource** | The authors' public browser repository describes knockout EOP, RB-TnSeq, phage genomes and BW25113/BL21 resources. This supports a separate source/provenance audit. | Inspect schemas, exact materials, assay meaning and publication status before selecting cases. Receptor dependence is not infection prediction or assembly evidence. |

PBIP source: [primary paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC12694490/), [pinned repository](https://github.com/a1678019300/PBIP/tree/a4c8a77032aaf7b474f1a316a350db207a6b508e). The inspected training scripts reconstruct scalers from training embeddings and train new networks; a checkpoint alone does not reproduce that pipeline. Its supplied strain split shares host IDs across training and test, so that split does not establish unseen-host performance. No checkpoint was loaded or target score calculated.

PhageMind source: [primary paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC13340165/), [pinned data description](https://github.com/YangSH-ac/PhageMind/blob/917b2ef3c62f8c730bf60b19fef1a896c6378b55/data/DATA.md). Archive-directory metadata includes BL21, `427_P2` and `T7_Portugal`, making overlap clearance a substantive requirement. Sequence and checkpoint payloads were not read in this triage.

The [pinned decoder](https://github.com/YangSH-ac/PhageMind/blob/917b2ef3c62f8c730bf60b19fef1a896c6378b55/code/src/my_model.py#L45) applies one linear output to the sum of host and phage embeddings. Algebraically, for fixed model and graph, `s(h,p) = a(h) + b(p) + c`; host-score differences therefore do not depend on the queried phage. This is an inference from the inspected implementation, not a numerical test or a claim that the paper fails. Query-dependent graphs or adaptation require separate analysis. The distinction matters when attributing host rankings to phage-dependent specificity.

Receptor-resource source: [author repository at inspected commit](https://github.com/mjohnson11/PhageDataSheets/tree/ae2b5e0e59b24f84666a8ebed01f5d1c0ded7b04/Ecoli_phages), associated with [the 2026 preprint](https://doi.org/10.64898/2026.04.02.716166). The README names `KEIO_EOP_reformatted.csv`, phage metadata and genetic-dependence data. Preprint status and incomplete full-text access remain explicit. These resources have not been admitted as training data or added to the frozen benchmark.

## Next concrete work

The [initial schema audit](../data/processed/receptor-resource-schema-audit.json) now confirms accessible files at the pinned author commit. The EOP table exposes `phage`, `Genotype`, `EOP` and `EOP (complemented)` fields; the metadata table provides genome accessions, separate BW25113/BL21 receptor and susceptibility fields, binding-protein annotations and source references. Only headers were inspected; full responses were retrieved, but outcome rows were not inspected or selected. These fields provide a concrete next source-audit path. Record-level provenance, experimental material identity, censoring rules and overlap remain to be established before admission.

1. Authenticate accessible receptor-resource metadata and schemas without choosing cases from favorable outcomes. Establish whether an additional mechanism panel can complement the original benchmark.
2. Audit the smallest complete PBIP inference bundle and its preprocessing/training membership. Resolve domain transfer before admitting it for native *E. coli* scoring; keep genus-only tools outside the strain-specific success criterion.
3. Select a method only after its input and mechanism requirements are supported. Freeze parameters, target exclusions, evaluation grouping and identical-subset baselines before producing scores.

An additional panel cannot replace the Ando question or turn receptor-class accuracy into strain-level or assembly validation. No new method is validated by this rerouting. The useful change is that project execution no longer waits on one inaccessible deposit.
