# An independent contrast source for method feasibility

Updated 2026-11-07. **The published BASEL–ECOR panel contains opposing host preferences that the two-host receptor-resource subset lacks.** It is a concrete candidate for secondary method validation, conditional on sequence, relatedness and endpoint admission. It does not replace the original Ando retrospective question.

## Why this source is being investigated

The [previous source-host audit](31-source-host-contrast-readiness.md) found two informative test comparisons sharing one host-preference direction. This permits a fixed host order to satisfy every comparison, making that subset insufficient to demonstrate a phage-dependent ranking advantage.

The previously authenticated [BASEL–ECOR roster](../data/processed/basel-validation-source-audit.json) names 52 phages, 25 ECOR hosts and 1,240 observed pairs. Its source is the [GenoPHI paper's experimental validation](https://doi.org/10.1038/s41564-026-02482-5), with a [pinned public matrix](https://github.com/Noonanav/GenoPHI/blob/f7e30abc1b0bbba1ae76b8f58a3363bbbc5cd54f/data/experimental_validation/BASEL_ECOR_interaction_matrix.csv). Phage accessions and host WGS master accessions were mapped before inspecting outcomes. These ECOR inputs are a different study panel from the unavailable Picard host-genome collection. No Picard archive, mirror or genome download was used for this audit.

## Executed outcome-distribution audit

The [policy](../software/tipcraft/benchmarks/receptor-resource/basel-ecor-contrast-policy.json) was fixed before this panel's first outcome-distribution calculation. The [script](../scripts/audit-basel-ecor-contrasts.py) retains **every observed pair and the complete 52/25 roster**. It rejects duplicate observations and values other than finite numeric 0/1. Unobserved combinations remain unknown. No source subset, score, classifier, calibration or target inference is selected or produced.

The [receipt](../data/processed/basel-ecor-contrast-audit.json) records:

| Source property | Observed result |
| --- | ---: |
| Phages / host identifiers | 52 / 25 |
| Explicit pair observations | 1,240 |
| Source 0 / 1 values | 938 / 302 |
| Missing combinations out of the 1,300 Cartesian pairs | 60 |
| Phage queries with both observed states | 49 |
| Within-query positive-negative host comparisons | 3,962 |
| Host pairs with at least one directional comparison | 291 |
| Host pairs with **both opposite directions**, across phages | **167** |

Each directional count uses distinct observed phages. For an unordered host pair, opposite directions mean that one phage has observed states `(1, 0)` while another has `(0, 1)` for those same host IDs. A single fixed host order cannot satisfy both. This establishes an opportunity to test interaction-dependent rankings; it does not show that any model succeeds.

The 3,962 comparisons are correlated, not independent experimental replicates. Phages, hosts and observation counts must be grouped before estimating uncertainty or testing generalization. The current audit has not applied genomic exclusions or a modelling split, so these numbers must not be quoted as final eligible test counts.

Bas01 and Bas02 have 25 observed zeros each. Bas14 has 24 observed zeros and one unknown combination. They remain in the panel rather than being removed for lacking a ranking contrast. They contribute no within-query positive-negative comparisons, but could matter under another properly specified endpoint.

## Admission boundaries

The source endpoint is **spot clearance versus no clearance**, not automatically productive infection or plaque formation. The previous source audit records the paper's exclusion of 60 unclear outcomes. This audit preserves those gaps; it does not fill them as zeros or independently certify a lower limit of detection. No mapping from these source values to quantitative EOP has been established.

The panel was already published as a validation resource for GenoPHI, and the method was selected with that literature available. Reusing it is exploratory retrospective component evaluation, not a fresh blind external test of the published method. The new contrast audit has now also exposed outcomes during development. Any later run must disclose that exposure and freeze its full method before generating scores.

Existing accession lists do not establish byte-identical assay-stock genomes or cleared homology. The published roster overlaps the original GenoPHI training hosts by name; any comparison to pretrained or paper-derived models must audit that overlap. A new independently partitioned fit would need its own complete manifest. BASEL IDs also occur in the receptor-resource collection; those records cannot be split into nominally independent training and testing phages by changing dataset names.

## Next concrete admission work

1. Resolve the 25 ECOR WGS master accessions to versioned assemblies and exact contig inventories from official sequence archives. Authenticate payloads and retain multi-contig genomes without substituting strain-name references.
2. Reconcile the 52 published phage accessions against authenticated source containers or official records. Do not assume that a shared BASEL name makes their sequences identical.
3. Apply predeclared target-related exclusions and grouped host/phage partitions before fitting. The current raw contrast counts must then be recomputed on the eligible evaluation subset without choosing favorable members.
4. Freeze a literature-backed method, endpoint, training-only feature construction, hyperparameters, baselines, attribution controls and uncertainty reporting. Compare a paired method with host-only and phage-only ablations on identical observed pairs; preserve unknown observations and zero-vector cases.

GO here means **continue source admission**, not train an admitted model or claim host prediction. Original Ando independent-method coverage remains **0/10**, and four exact exchange-construct cases still abstain. Binding, assembly and functional infection remain separate questions.

## Other source leads and rejected shortcuts

The [Viral Host Range database source 135](https://viralhostrangedb.pasteur.cloud/data-source/135/) provides a separate Mt1B1/ECOR spotting-data lead. Its descriptive text says 17 phages/90 strains, while its current inventory lists 28/98 and three named phage accessions. That discrepancy, the downloadable payload and the primary publication correspondence require checking before admission. The page's generic T7/example matrix is **not** imported as a study observation.

The [BioML hackathon deposit](https://zenodo.org/records/13956128) describes its E. coli archive as derived from the same Picard/Guelin study. It is therefore not a new independent source and is not downloaded as a workaround for the skipped host collection. The existence of a different download URL does not create independent evidence.

## Reproduction

```bash
uv run --project software/tipcraft python scripts/audit-basel-ecor-contrasts.py --output /tmp/basel-ecor-contrast-audit.json
```

Four tests cover opposing directions, missing observations, uninformative all-positive queries, duplicate observations, nonbinary values and nonfinite values. The complete 279-test suite passes with Biopython 1.86 enabled.
