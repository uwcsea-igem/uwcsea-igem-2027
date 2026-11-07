# Official BASEL accessions and author-cache correspondence

Updated 2026-11-07. **All 52 official phage records have been retrieved and verified across annotation/FASTA formats. Forty-nine exactly match the existing author-cache DNA; three do not.** Official-record authentication and correspondence to an older cache are separate claims. The three mismatches require new relatedness work before their official genomes can inherit any cached group assignment.

## Frozen full-roster execution

The [policy](../software/tipcraft/benchmarks/receptor-resource/basel-phage-accession-policy.json) was recorded before bulk retrieval and binds the prior 52-ID/accession authority, authenticated author-container inventory and [three-record route pilot](34-ecor-host-groups-and-target-exclusions.md#bounded-phage-record-route-pilot). The pilot bytes for Bas01, Bas12 and Bas14 are explicitly reused; the other 49 accessions are retrieved under the same prefix-specific rules.

The [producer](../scripts/authenticate-basel-phage-accessions.py) uses ENA EMBL/FASTA for the 41 PZ and nine MZ accessions, and NCBI EFetch GenBank/FASTA for the two RefSeq NC accessions. Biopython 1.86 must parse exactly one nonempty IUPAC DNA record per official format, with an explicit versioned accession matching the request. The two official formats must agree on both accession version and exact upper-case DNA. Original author GenBank file hashes are authenticated before comparison; custom record names are not accepted as sequence-identity proof.

Every source ID is retained, including failed correspondences. No phenotype table is loaded. GenBank/EMBL annotation records are parsed to retrieve accession identity and DNA; receptor annotations are not used for selection or matching. No model, score or outcome-selected subset is produced. There is no sequence trimming, concatenation, rotation/reverse-complement search, editing or nearest-reference substitution.

## Two distinct results

The original [execution receipt](../data/processed/basel-phage-accession-authentication.json) applies the strict combined official-plus-author-correspondence gate. It therefore records **49 authenticated exact correspondences and three unresolved correspondences**. That receipt remains unchanged.

A separate [offline replay](../scripts/replay-basel-phage-accessions.py) re-parses all saved official payloads and author records and makes the two claims explicit. Its [receipt](../data/processed/basel-phage-accession-replay.json) establishes:

| Claim | Verified scope |
| --- | ---: |
| Complete official annotation/FASTA pairs | 52 / 52 |
| Official records resolving to sequence version `.1` | 52 / 52 |
| Total bases across the official genomes | 4,950,501 |
| Exact official-versus-author DNA matches | 49 / 52 |
| DNA correspondences not exact | 3 / 52 |
| Raw artifacts re-authenticated in replay | 156 |

All official format pairs agree, including the three records whose author-cache correspondence fails. The difference is therefore not a missing official record, a failed network request or a disagreement between the two downloaded official formats.

| Phage | Official version | Official length | Author-cache length |
| --- | --- | ---: | ---: |
| Bas06 | PZ445410.1 | 51,474 | 51,473 |
| Bas21 | PZ445422.1 | 56,589 | 56,738 |
| Bas68 | PZ445457.1 | 39,466 | 39,465 |

Both versions and their hashes are retained. These length differences are descriptive input facts, not an alignment, inferred edit history, sequencing-error diagnosis or assessment of biological function. No attempt was made to repair either representation or infer which experimental stock it matches. Equality of the other 49 DNA strings does not establish equality of protein annotations or physical assay stocks either.

## Consequence for genomic grouping and modelling

The previous [whole-phage groups](27-phage-groups-and-source-split.md) describe the authenticated author-container DNA. For the 49 exact matches, that exact sequence correspondence supports reuse of the existing **genotype** evidence, subject to the existing exclusions and grouping rules. It is not training admission or proof of an unseen group.

The official Bas06, Bas21 and Bas68 records must be treated as additional sequence representations. They cannot silently inherit their namesakes' numeric similarities or component IDs. A subsequent frozen comparison must retain original targets and relevant source nodes as bridges, apply the same project thresholds and account for both official and author versions explicitly. The old matrix, mismatch receipts and original case observations remain preserved.

The next step is to resolve these three additional genotypes against the existing reference graph and then compute the eligible secondary-panel manifest and contrasts. A label-free sequence mismatch is not grounds to drop a difficult case for better performance. If a version cannot be admitted, it retains an explicit reason and denominator; if official-reference analysis proceeds, its new grouping must be measured. No smaller subset is presented here as the complete validated panel.

The [23 candidate host references in 15 groups](34-ecor-host-groups-and-target-exclusions.md) remain unchanged. Complete source endpoint, stock/reference, partition and pre-score method requirements still apply. Original Ando independent-method coverage remains **0/10**; four exact exchange-construct cases still abstain. No binding, assembly or infection success is inferred from this source authentication.

## Verification and reproduction

The full replay verifies 104 official payloads and 52 author containers, reproduces every exact-match/mismatch status and distinguishes official version identity from author-cache identity without any network request or model execution. Four tests cover exact DNA equality, version disagreement, ENA header resolution, missing/wrong accessions, multiple records and invalid DNA. The complete 293-test suite passes with Biopython 1.86 enabled.

With the authenticated source/pilot caches and fresh outputs:

```bash
uv run --project software/tipcraft --with biopython==1.86 python scripts/authenticate-basel-phage-accessions.py --workdir data/raw/genophi-external-validation/basel-phage-accessions-reproduction-1 --output /tmp/basel-phage-accession-authentication.json
uv run --project software/tipcraft --with biopython==1.86 python scripts/replay-basel-phage-accessions.py --report /tmp/basel-phage-accession-authentication.json --output /tmp/basel-phage-accession-replay.json
```

Retrieval timestamps and local paths may differ. The fixed source hashes, explicit accession versions and per-sequence hashes are the reproducibility controls. Future upstream changes must not silently overwrite the executed input set.
