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
| **Camejo et al. 2025 PPI-based models** | The primary study explicitly addresses strain-specific outcomes using inferred domain-interaction features. | Models are fitted separately per phage, using growth-inhibition labels. The inspected work does not establish arbitrary unseen-phage transfer to T7/T3 or equivalence to plaque formation. |

PBIP source: [primary paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC12694490/), [pinned repository](https://github.com/a1678019300/PBIP/tree/a4c8a77032aaf7b474f1a316a350db207a6b508e). The inspected training scripts reconstruct scalers from training embeddings and train new networks; a checkpoint alone does not reproduce that pipeline. Its supplied strain split shares host IDs across training and test, so that split does not establish unseen-host performance. No checkpoint was loaded or target score calculated.

PhageMind source: [primary paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC13340165/), [pinned data description](https://github.com/YangSH-ac/PhageMind/blob/917b2ef3c62f8c730bf60b19fef1a896c6378b55/data/DATA.md). Archive-directory metadata includes BL21, `427_P2` and `T7_Portugal`, making overlap clearance a substantive requirement. Sequence and checkpoint payloads were not read in this triage.

The [pinned decoder](https://github.com/YangSH-ac/PhageMind/blob/917b2ef3c62f8c730bf60b19fef1a896c6378b55/code/src/my_model.py#L45) applies one linear output to the sum of host and phage embeddings. Algebraically, for fixed model and graph, `s(h,p) = a(h) + b(p) + c`; host-score differences therefore do not depend on the queried phage. This is an inference from the inspected implementation, not a numerical test or a claim that the paper fails. Query-dependent graphs or adaptation require separate analysis. The distinction matters when attributing host rankings to phage-dependent specificity.

Receptor-resource source: [author repository at inspected commit](https://github.com/mjohnson11/PhageDataSheets/tree/ae2b5e0e59b24f84666a8ebed01f5d1c0ded7b04/Ecoli_phages), associated with [the 2026 preprint](https://doi.org/10.64898/2026.04.02.716166). The README names `KEIO_EOP_reformatted.csv`, phage metadata and genetic-dependence data. The resource is a preprint; its v1 Methods and supplementary figures have now been inspected. These resources have not been admitted as training data or added to the frozen benchmark.

## Paired-host encoder sanity check

The [two-reference host encoder test](26-host-encoder-capacity.md) now provides an executed negative applicability result before interaction-model fitting: all three original target-host proteomes map to the same nonzero vector under the proposed training-pattern-merged5-mer union projection. This candidate cannot act as the sole host representation in a sequence-feature paired ranker for this target panel. This does not reject the receptor-only endpoint or all GenoPHI implementations, and no new host-range accuracy is claimed.

## Additional strain-method applicability check

[Camejo et al. 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12579225/) construct protein-domain interaction features and fit separate models for individual phages. Their endpoint is growth inhibition; it is not automatically a plaque label or direct receptor binding. The authors also discuss sensitivity to phage mutations and the possible need for phage-specific models. Applying this recipe by fitting on the Ando T7/T3 outcomes would violate this project's target-label exclusion. An independently justified unseen-phage adaptation and authenticated training/inference bundle would be needed before admission. The inspected availability statement lists deposited phage genomes; it does not authenticate a runnable T7/T3 predictor. This is a bounded applicability finding, not a conclusion about every external asset or the paper's performance.

## Next concrete work

### Phage-only receptor-annotation modelling

The receptor study supplies a method route that does not require the Picard host genomes. Its [v1 Methods, Host target predictive modelling](https://www.biorxiv.org/content/10.64898/2026.04.02.716166v1.full) describes GenoPHI models using amino-acid 5-mers from all phage proteins and 13 separate binary receptor-dependency tasks: Tsx, OmpA, OmpF, OmpC, FhuA, BtuB, LptD, LamB; Kdo, HepI, HepII, GluI; and NGR (N4-glycan receptor polysaccharide). A dependency label is not automatically direct-binding evidence.

The source model-role column contains 206 training, 49 testing and five mutant entries. T7 is explicitly in training. The paper's final models incorporate the 49 prospective cases, so those cases cannot be called unseen when evaluating those final models. A new retrospective fit must freeze its own sequence-related exclusions and evaluation membership. The metadata-only roster uniquely maps all 260 primary IDs to GenBank filenames in a 313-file directory. These candidate payloads still require authentication and a sequence manifest; the 53 catalogue extras are not automatically admitted.

The pinned GenoPHI [de novo workflow](https://github.com/Noonanav/GenoPHI/blob/f7e30abc1b0bbba1ae76b8f58a3363bbbc5cd54f/genophi/workflows/kmer_full_workflow.py#L510) ignores protein families; its [k-mer table implementation](https://github.com/Noonanav/GenoPHI/blob/f7e30abc1b0bbba1ae76b8f58a3363bbbc5cd54f/genophi/workflows/kmer_table_workflow.py) builds binary occurrence features and merges matching patterns. The existing novel-genome assignment workflow instead uses a protein-family join. Therefore a direct global 5-mer projection needs an explicit, independently verified adaptation. Exact paper commands, seeds and split manifests have not been established.

This is a **conditional route for an additional receptor-annotation backend**. Its output alone cannot rank BL21, BW25113 and MG1655: that would additionally require a justified host-context model and independent validation. The original strain-ranking objective and assembly gate remain separate. No model has been fitted or admitted by this source review. The [readiness assessment](24-receptor-method-readiness.md) now records source-label version differences, direct 5-mer implementation boundaries and the remaining admission gates.

The [record-level source audit](23-receptor-resource-provenance.md) now authenticates the two tables and the official v1 supplementary workbook. It retains 923 browser records and now verifies all 1,846 original/complemented values against Dataset S4 within a preregistered serialization tolerance. Official figures establish plotted log10 EOP; per-cell scale/censoring and quantitative training admission remain separate gates.

The [initial schema audit](../data/processed/receptor-resource-schema-audit.json) now confirms accessible files at the pinned author commit. The EOP table exposes `phage`, `Genotype`, `EOP` and `EOP (complemented)` fields; the metadata table provides genome accessions, separate BW25113/BL21 receptor and susceptibility fields, binding-protein annotations and source references. That initial report inspected headers only; subsequent row-preserving audits are documented separately above. These fields provide a concrete next source-audit path. Record-level provenance, experimental material identity, censoring rules and overlap remain to be established before admission.

1. Complete container-aware whole-genome grouping with frozen exclusions and original-R numeric validation. Keep source-label versions, annotation conventions and all unknown states explicit.
2. Authenticate the newly located author host-genome/protein assets and establish what stocks/references they represent. Assess a paired host-feature adaptation separately from a phage-only receptor classifier; two source host backgrounds do not establish general unseen-host validity. Production hosts never fill missing MG1655 labels.
3. Select a scoring method only after its input, endpoint and mechanism requirements are supported. Freeze parameters, target exclusions, evaluation grouping and identical-subset baselines before producing scores. An applicable model may be tested and fail; a source or task mismatch remains HOLD rather than invented coverage.
4. Retain PBIP as a conditional asset-audit alternative. Resolve authenticated preprocessing/training membership and *Klebsiella*-to-*E. coli* transfer before admitting it for native strain scoring.

An additional panel cannot replace the Ando question or turn receptor-class accuracy into strain-level or assembly validation. No new method is validated by this rerouting. The useful change is that project execution no longer waits on one inaccessible deposit.
