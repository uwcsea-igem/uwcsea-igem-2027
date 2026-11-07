# Phage relatedness groups and independent holdout boundaries

Updated 2026-11-07. **The alternative phage collection now has a complete, original-R-checked relatedness matrix. Its published training/testing roles overlap under the project's genomic grouping rule.** A deterministic candidate split preserves all 49 eligible testing phages and quarantines their related training phages; it does not yet admit labels, features or a model.

## Completed alignment and numeric validation

The unchanged BLAST search processed **270 separate nucleotide records**: 268 records from 261 source containers, plus two original target-reference controls. Its raw output contains 1,245,085 HSPs. The original v1 receipt preserves the parser's failure on legitimate multiple completion batches. A separately frozen [protocol-only replay](../data/processed/receptor-whole-genome-alignments-v2.json) validates all queries in six batches (`55, 56, 46, 47, 43, 23`) without rerunning BLAST or changing HSP order, strings, coordinates or numeric literals.

The [v2 grouping policy](../software/tipcraft/benchmarks/receptor-resource/whole-genome-group-policy-v2.json) retains the original scientific rules. Resource-only changes stream the external calculator's input and increase predeclared process limits before numeric execution. The [allocation control](../data/processed/viridic-streaming-allocation-control.json) reproduces the previously R-verified 99-genome result exactly, including all 4,851 pair records and 99 self-controls. The upstream numeric functions remain unchanged.

For multiple-record containers, searches remain separate. Cumulative offsets transform coordinates solely for numeric overlap accounting; no DNA is concatenated, trimmed or invented. Every HSP stays within its original record. The [executed grouping report](../data/processed/receptor-whole-genome-groups.json) and [full matrix](../data/processed/receptor-whole-genome-similarities.tsv) establish:

| Scope | Executed result |
| --- | ---: |
| Whole-container nodes | 263 |
| Complete unordered pairs | 34,453 |
| Nonself directions compared with original R | 68,906 |
| Maximum unrounded difference (percentage points) | 5.684341886080802e-14 |
| Directional counts, span diagnostics and three-decimal scores | Exact agreement |
| Relatedness components | 52 |

All 263 dataset-specific self-controls pass the rule frozen before calculation. The generic helper's all-base self-check **still records a failure for TP15**: both implementations report 87,528 self-identical bases for its 87,728-base reference, which contains 200 Ns. The separate, predeclared ambiguous-sequence rule requires exact R/Python agreement and finite bounded counts rather than forced 100% self-identity. The generic failure is preserved in the report; it was not silently changed into a generic pass. This establishes numeric agreement on this dataset, not universal equivalence or biological stock identity.

The CPU-only search took about 68.8 minutes; external Python calculation about 102 seconds and original R about 376 seconds. Hash verification and transformation add separate work. Reported resource counters are serial direct-child measurements and cumulative child peak RSS, not aggregate workflow memory.

## Target-related exclusions

Unrounded similarity of at least **69.99%** forms an undirected edge: the predeclared 70% threshold minus the 0.01-point exclusion buffer. All source, mutant, catalogue-extra and original-reference nodes remain as possible bridges. These components are conservative project relatedness units, not taxonomic assignments.

| Target-containing component | Source members excluded |
| --- | --- |
| Original T7 reference | Bas102, Bas103, Bas64, Bas65, T7 |
| Original T3 reference | Bas66, catalogue-extra T3 |

Six excluded source members belong to the published training cohort; catalogue T3 is outside the cohort. Direct gp17 near-hit seeds add no further components in this dataset. Bas64 is excluded through whole-genome connectivity, independently of its separately recorded broad gp17 flag. The catalogue T3 remains distinct from the original target accession and is never substituted for an Ando construct.

After target exclusions, 254 cohort records remain in 50 components: 200 published training members, 49 testing members and five mutants. These are sequence candidates, not admitted training observations.

## Published roles versus genomic independence

The [source-role audit policy](../software/tipcraft/benchmarks/receptor-resource/source-role-overlap-policy.json) fixes a label-free rule after grouping and before any model fit. The [executed audit](../data/processed/receptor-source-role-overlap.json) finds **nine components shared by the eligible published training and testing sets**. They contain 42 of the 49 testing phages and 116 of the 200 training phages.

Prospectively isolated or experimentally assayed phages can still be close genomic relatives of training phages. This finding concerns the project's grouping criterion; it does not dispute the source study's stated prospective-isolate validation. This project cannot describe that original split as genomic-group-disjoint.

To retain every eligible source test member, the candidate partition quarantines all training members in any testing-containing component:

| Cohort accounting | Records | Distinct fitting/evaluation groups |
| --- | ---: | ---: |
| Target-related exclusions | 6 | Not fitting/evaluation candidates |
| Candidate training | 84 | 34 |
| Candidate testing, all retained | 49 | 16 |
| Related training members quarantined | 116 | Not fitted |
| Mutants, nontraining bridge records | 5 | Not fitted |

All 260 cohort records are accounted for exactly once. Candidate training and testing share zero genomic groups. No random fallback, phenotype-based selection or metric-driven removal occurs. Source labels were already visible during method development: this remains exploratory retrospective work, not a newly blinded test. Endpoint-specific label availability, rare classes, source versions, host contexts and feature admission remain unresolved gates.

## Decision consequence

The grouping evidence supports a reproducible independent-validation boundary. It does not supply a host-range predictor. The [host-encoder collision](26-host-encoder-capacity.md) still makes the proposed merged/union host representation NO-GO for the current paired ranker, including under uniform annotation. Adding phage samples or changing the split cannot repair identical host vectors while keeping that host encoder fixed. Additional host variation or a different representation would require a new fit/check protocol.

No model has been fitted on this candidate partition. Original Ando independent-method coverage remains **0/10**; the WIsH failure and four exact-construct abstentions remain intact. Any later method needs an explicit pre-score contract and matched-subset baselines. Receptor annotation, binding, assembly and functional infection remain separate endpoints.

## Reproduction

With authenticated caches and the pinned runtimes, use fresh outputs:

```bash
uv run --project software/tipcraft python scripts/aggregate-receptor-genomes.py --alignment-report data/processed/receptor-whole-genome-alignments-v2.json --workdir data/raw/receptor-resource/whole-genome-groups-reproduction-1 --output /tmp/receptor-whole-genome-groups.json --matrix /tmp/receptor-whole-genome-similarities.tsv
uv run --project software/tipcraft python scripts/audit-receptor-role-overlap.py --output /tmp/receptor-source-role-overlap.json
```

The role audit consumes the versioned verified group/matrix artifacts; it does not rerun alignments or read phenotype labels. The preserved v1/v2 execution lineage and the original-R generic/dataset-specific results remain inspectable separately.
