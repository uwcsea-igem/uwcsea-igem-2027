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

The [official v1 supplementary workbook](https://doi.org/10.6084/m9.figshare.31930314.v1) was actually downloaded and checked against its published size and MD5. Dataset S4 has 1,846 literal measurement rows: **1,404 negative, 17 zero and 425 positive**. Its header labels the field EOP without establishing a transformation. These totals agree with the two browser contexts combined; exact row-wise correspondence has not yet been verified. The [source follow-up report](../data/processed/receptor-measurement-source-followup.json) binds the workbook and observations. No logarithm or inverse transformation is applied.

The [author formatter](https://github.com/mjohnson11/PhageDataSheets/blob/ae2b5e0e59b24f84666a8ebed01f5d1c0ded7b04/Ecoli_phages/DataAssembly/format_EOP.ipynb) copies an external CSV's values unchanged and inner-merges the two contexts. It does not document the earlier value generation; unpaired source contexts could be omitted. The browser's colour range and numeric display do not establish units, detection limits or the meaning of missing annotations. This finding concerns interpretation of the stored data, not a claim that the experiment is wrong.

## Reproduction and next admission gate

```bash
uv run --project software/tipcraft python scripts/audit-receptor-resource.py --cache data/raw/receptor-resource --workdir data/raw/receptor-resource/reproduction-1 --output /tmp/receptor-resource-audit.json
```

Use fresh output locations and the authenticated cached files. Row-level raw values and source locators remain in ignored storage; the versioned report binds them by hashes. Tests cover nonfinite and negative entries, explicit bounds, missing cells, duplicate/ambiguous identities, schema errors and source tampering.

Before quantitative use, resolve the stored scale and censoring, verify source/material correspondences and freeze each case's role and relatedness exclusions. A mechanism panel must preserve empty-vector and complemented contexts and remain separate from the original Ando evaluation. No model has been fitted, no new target score exists and neither binding nor assembly is inferred.
