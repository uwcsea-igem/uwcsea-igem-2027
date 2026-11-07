# Receptor-resource provenance and measurement context

Updated 2026-11-07. **The alternative route now has authenticated public tables and a reproducible, row-preserving audit.** The audit connects source identifiers while separating literal values, experimental context and biological interpretation. These resources are independent of the declined Picard whole-host download.

## Authenticated source inventory

The [source policy](../software/tipcraft/benchmarks/receptor-resource/source-audit-policy.json) was written before outcome parsing. It fixes both files at author commit `ae2b5e0e59b24f84666a8ebed01f5d1c0ded7b04`, requires exact byte sizes/Git blob hashes and retains every row. No favorable-outcome selection, averaging or binary infection conversion occurs.

The [executed audit](../data/processed/receptor-resource-audit.json) records:

| Scope | Result |
| --- | ---: |
| Measurement rows / literal phage IDs / genotype literals | 923 / 140 / 38 |
| Exact primary-ID links to metadata | 909 rows |
| Unlinked rows | 14: four `CPL00121`, ten `MPec26` |
| Duplicate phage/genotype pairs or schema-issue rows | 0 |
| Metadata records | 260 |
| Nonblank accession / reference fields | 260 / 219 |

Identifier links and populated accession fields do not authenticate stocks or externally validate accession mappings. Six measurement rows have an exact T7 name warning. The absence of an exact T3 name is not sequence-homology clearance. Receptor and susceptibility annotations are phenotype-informed source metadata, not independent predictions.

## Literal values and source interpretation

| Browser field | Nonnegative numeric literals | Negative literals preserved as uninterpreted | Blank / explicit inequality cells |
| --- | ---: | ---: | ---: |
| EOP | 116 | 807 | 0 / 0 |
| EOP (complemented) | 326 | 597 | 0 / 0 |

Every row remains unadmitted, including nonnegative or zero entries. These are parser categories, not valid/invalid experiment counts or infection classes. Zero cannot be assigned a biological meaning until the stored scale is established.

The [v1 preprint](https://www.biorxiv.org/content/10.64898/2026.04.02.716166v1.full), Methods: Phenotypic validations, defines EOP as a mutant-to-parental plaque-count ratio. The contexts involve BW25113-derived KEIO deletions with empty-vector controls and corresponding ASKA complementation. Complementation is a separate observation, not assumed rescue. Methods describe technical duplicates, whereas the Dataset S4 legend describes biological duplicates; the source discrepancy remains recorded.

The [official v1 supplementary workbook](https://doi.org/10.6084/m9.figshare.31930314.v1) was actually downloaded and checked against its published size and MD5. Dataset S4 has 1,846 literal measurement rows: **1,404 negative, 17 zero and 425 positive**. Its header labels the field EOP without establishing a transformation. These totals agree with the two browser contexts combined. The subsequent frozen correspondence audit below verifies every source row against its browser context. The [source follow-up report](../data/processed/receptor-measurement-source-followup.json) binds the workbook and observations. No logarithm or inverse transformation is applied.

The [author formatter](https://github.com/mjohnson11/PhageDataSheets/blob/ae2b5e0e59b24f84666a8ebed01f5d1c0ded7b04/Ecoli_phages/DataAssembly/format_EOP.ipynb) copies an external CSV's values unchanged and inner-merges the two contexts. It does not document the earlier value generation; its inner merge could omit unpaired contexts in general. The executed comparison finds no such omission relative to the published Dataset S4. The browser's colour range and numeric display do not establish units, detection limits or the meaning of missing annotations. This finding concerns interpretation of the stored data, not a claim that the experiment is wrong.

## Executed correspondence and scale evidence

The [correspondence policy](../software/tipcraft/benchmarks/receptor-resource/measurement-correspondence-policy.json) fixed the exact keys and absolute serialization tolerance `1e-8` before comparison. The [executed report](../data/processed/receptor-measurement-correspondence.json) independently parses the workbook's OOXML and matches all **1,846 rows** to 923 original and 923 complemented browser contexts. Seventeen numeric pairs are Decimal-exact; 1,829 differ within the fixed tolerance. There are zero mismatches, ambiguous keys, unrecognized suffixes or unpaired rows. The maximum absolute difference is `4.822699E-9`. No numerical transform or phenotype classification is involved.

The official v1 **Figures S2 and S3 explicitly label their EOP axes “Efficiency of plating (log10)”**. Both downloaded images were verified against the deposit's size and MD5 and visually inspected; the [figure evidence record](../data/processed/receptor-measurement-scale-evidence.json) binds their SHA256 hashes. This establishes the plotted scale and supports interpreting the associated negative values as logarithmic measurements, rather than absence of infection. Row correspondence establishes the workbook/browser serialization bridge. It does not establish a universal detection floor or independently certify the units and censoring of every tabulated cell. Those remain required for quantitative admission.

The initial literal audit and its uninterpreted-negative categories remain unchanged as the record of what was established at that stage. Neither positive nor negative log10 EOP alone is converted to a binary infection label. The technical/biological replicate discrepancy and unlinked identifiers remain open; `CPL00221` shown in Figure S3 is not silently equated with `CPL00121` in the table.

The [metadata-only genome roster](../data/processed/receptor-phage-roster-audit.json) uniquely matches all 260 primary IDs to pinned GenBank filenames: 206 training, 49 testing and five mutant entries. The directory contains 313 GenBank files, including 53 outside that cohort. Filename correspondence is not sequence or stock authentication. T7 is in the source training cohort; T3.gbk is an extra catalogue file. Payload authentication and sequence-related exclusion precede any fit.

## Existing T7/T3 container authentication

The [frozen two-container policy](../software/tipcraft/benchmarks/receptor-resource/t7-t3-container-policy.json) selects the original target names for input inspection, without receptor-outcome filtering. The [executed container audit](../data/processed/receptor-target-container-audit.json) authenticates both files against their pinned Git blob hashes and parses existing annotations with Biopython 1.86. No sequences or translations were designed or filled in.

| Catalogue file | Deposited record | Length | CDS with deposited translation | CDS missing translation |
| --- | --- | ---: | ---: | ---: |
| T7.gbk | NC_001604.1 | 39,937 | 61 | 0 |
| T3.gbk | KC960671.1 | 38,209 | 54 | 0 |

Both parse without warnings and declare translation table 11 for their CDS features. Their deposited translations contain standard amino-acid letters; translation agreement with the nucleotide annotations has not yet been checked. The T7 sequence hash exactly matches the [original WIsH reference manifest](../software/tipcraft/benchmarks/retro-t7-t3/wish-inputs.json). **The catalogue T3 is a different accession and sequence hash from the original NC_003298.1 reference (38,208 bp).** It cannot silently replace that target input or establish identity with the Ando experimental stock. This audit therefore advances input feasibility without increasing the original benchmark's prediction coverage.

```bash
uv run --project software/tipcraft --with biopython==1.86 python scripts/audit-receptor-target-containers.py --cache data/raw/receptor-resource/target-containers-v1 --output /tmp/receptor-target-container-audit.json
```

The separate [method-readiness assessment](24-receptor-method-readiness.md) compares the official v1 metadata labels with the updated browser table. Their 260 IDs and model roles agree, but 12 receptor/LPS annotation cells differ across 11 phages. Source version therefore remains part of the model specification.

## Reproduction and next admission gate

```bash
uv run --project software/tipcraft python scripts/audit-receptor-resource.py --cache data/raw/receptor-resource --workdir data/raw/receptor-resource/reproduction-1 --output /tmp/receptor-resource-audit.json
```

The separate correspondence audit is reproduced with:

```bash
uv run --project software/tipcraft python scripts/compare-receptor-measurements.py --workdir data/raw/receptor-resource/correspondence-reproduction-1 --output /tmp/receptor-measurement-correspondence.json
```

Use fresh output locations and the authenticated cached files. Row-level raw values and source locators remain in ignored storage; the versioned report binds them by hashes. Tests cover nonfinite and negative entries, explicit bounds, missing cells, duplicate/ambiguous identities, schema errors and source tampering.

Before quantitative use, resolve the stored scale and censoring, verify source/material correspondences and freeze each case's role and relatedness exclusions. A mechanism panel must preserve empty-vector and complemented contexts and remain separate from the original Ando evaluation. No model has been fitted, no new target score exists and neither binding nor assembly is inferred.
