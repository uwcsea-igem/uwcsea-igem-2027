# Versioned ECOR reference inputs

Updated 2026-11-07. **All 25 ECOR host inputs in the published secondary panel now have authenticated official WGS contig sets.** This resolves a reference-input prerequisite independently of the skipped Picard host collection. It does not establish assay-stock equivalence, genetic independence or predictive performance.

## Roster and route

The [prior source audit](../data/processed/basel-validation-source-audit.json) maps every observed-panel host to its WGS master accession before outcome-distribution inspection. That complete 25-host roster is unchanged. The supplementary assembly table also contains ECOR36 and ECOR37, but neither is an observed-panel host; they were not added or selected using phenotypes.

Official [ENA retrieval documentation](https://ena-docs.readthedocs.io/en/latest/retrieval/programmatic-access/browser-api.html) and [contig-set documentation](https://ena-browser-docs.readthedocs.io/en/latest/browser/ena-record.html) define the master-record and complete-set download routes. A route pilot inspected ECOR12's master record, assembly XML and compressed FASTA before the bulk policy was recorded. The ECOR12 assembly XML identifies `GCA_003334775.1`, `SAMN04158349`, project `PRJNA230969` and its WGS FASTA link. That independent XML cross-check is specific to the pilot; GCA accessions were not inferred for the other 24 hosts.

The [bulk policy](../software/tipcraft/benchmarks/receptor-resource/ecor-wgs-input-policy.json) discloses the pilot and reuses its hash-authenticated master/FASTA bytes. The remaining 24 hosts use the same official route. The [retriever](../scripts/authenticate-ecor-wgs-inputs.py) reads only the accession authority and assembly metadata, not interaction outcomes. It retains every roster entry with an explicit status, rejects malformed or incomplete sets and records retrieval locations and hashes.

## Executed result

The [authentication receipt](../data/processed/ecor-wgs-input-authentication.json) establishes:

| Scope | Result |
| --- | ---: |
| Frozen roster / authenticated hosts | 25 / 25 |
| Unresolved inputs | 0 |
| Complete contig records retained | 3,477 |
| Total nucleotide bases | 120,989,458 |
| Ambiguous nucleotide characters | 0 |
| Master counts disagreeing with source `Scaffolds` | 0 |
| Compressed FASTA payload bytes, including reused pilot | 37,911,609 |

Every master belongs to project `PRJNA230969` and declares its expected `MOD1-ECORn` strain. All 25 WGS builds resolve to `01`, with 25 distinct BioSample accessions. All downloaded contig records have sequence version `.1`. **WGS build and contig sequence version are separate identifiers**; neither is inferred from the other.

Each input must contain exactly the master's complete inclusive contig range, without duplicates or unexpected accessions. The expanded payload must pass gzip integrity checking and contain nonempty IUPAC DNA with matching versioned accession and strain declarations. Every contig's accession, sequence version, length and sequence SHA-256 is retained in a full raw manifest. No sequence is concatenated, trimmed, filled or reduced to a longest contig. The public receipt binds every raw master, compressed/expanded FASTA and full manifest.

## Offline replay and independence limits

The [offline replay](../scripts/replay-ecor-wgs-inputs.py) verifies all **100 raw artifact hashes**, re-expands each gzip, re-parses every master and every contig, and compares the complete inventories and public counts. Its [receipt](../data/processed/ecor-wgs-input-replay.json) records the same 25 hosts, 3,477 contigs and 120,989,458 bases without a network request, gene caller or model execution.

A separately defined exact fingerprint hashes the sorted multiset of `(contig length, sequence SHA-256)` tuples, preserving repeated contig multiplicity. No two inputs have identical fingerprints under this rule. This check is not reverse-complement or fragmentation invariant, and it does **not** establish that all 25 references are genetically independent. Near-relatedness grouping remains mandatory.

The matching source assembly metrics are corroborating metadata, not experimental-stock authentication. These are public reference representations tied to the paper's accession list. A strain name, BioSample or assembly agreement does not prove that the reference is byte-identical to the host stock used in a later susceptibility assay. No susceptibility annotation has been reassigned to a substitute genome.

## Consequence and next step

GO for **reference-input preparation and relatedness checks**. The next admission work is to group these hosts with the original target references under the previously specified sequence-overlap rules, and reconcile the 52 phage accessions against exact public sequence records. Complete an eligible, group-aware evaluation manifest before fitting or scoring.

The subsequent [complete host comparison](34-ecor-host-groups-and-target-exclusions.md) now supplies that reference-group boundary: ECOR13 and ECOR25 are excluded as target-related, leaving 23 candidate references in 15 groups. This does not turn input authentication or grouping into final training admission.

The [raw panel's opposing host preferences](32-secondary-contrast-source-selection.md) remain a source-feasibility finding. This input authentication adds no prediction, comparison against a baseline, binding evidence or assembly evidence. Original Ando independent-method coverage remains **0/10**; the measured WIsH failure and four exact-construct abstentions are preserved.

## Reproduction

The bulk policy binds two ECOR12 route-pilot cache files. On a fresh checkout, recover those exact official bytes first; hashes are checked before use:

```bash
mkdir -p data/raw/genophi-external-validation/ena-metadata-probe-v1
curl --fail --location https://www.ebi.ac.uk/ena/browser/api/embl/QOWX00000000 --output data/raw/genophi-external-validation/ena-metadata-probe-v1/master.txt
curl --fail --location https://ftp.ebi.ac.uk/pub/databases/ena/wgs/public/qow/QOWX01.fasta.gz --output data/raw/genophi-external-validation/ena-metadata-probe-v1/QOWX01.fasta.gz
uv run --project software/tipcraft python scripts/authenticate-ecor-wgs-inputs.py --workdir data/raw/genophi-external-validation/ecor-wgs-inputs-reproduction-1 --output /tmp/ecor-wgs-input-authentication.json
uv run --project software/tipcraft python scripts/replay-ecor-wgs-inputs.py --report /tmp/ecor-wgs-input-authentication.json --output /tmp/ecor-wgs-input-replay.json
```

Use fresh output destinations and the authenticated source accession/assembly tables. Updated upstream bytes must fail the frozen hash check rather than silently replacing the pilot. Retrieval timestamps and local paths can differ between executions; the fixed input identities and per-sequence hashes are the reproducibility controls.

Six focused tests cover complete-set membership, per-record sequence versions, ambiguity preservation, partial/duplicate/unexpected records, incorrect master/project/strain declarations and invalid DNA. The complete 285-test suite passes with Biopython 1.86 enabled.
