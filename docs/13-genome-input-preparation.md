# Genome inputs and relatedness audit

Updated 2026-11-07. **Input preparation has advanced from file-name candidates to authenticated public genome downloads and cross-source label curation.** The next prediction experiment remains unrun while complete host inputs, phage exclusions and the training specification are prepared.

## Current evidence

| Input task | Measured result | Consequence |
| --- | --- | --- |
| Original-to-GenoPHI label reconciliation | 37,632 nonblank pairs agree; 156 original-source-empty pairs appear downstream as zero | Exclude all 156 from our future run pending clarification. [Raw-annotation review](14-training-label-provenance.md) establishes ambiguity, not replacement labels. |
| Host naming correspondence | All 402 model host IDs have a unique punctuation-normalized original-table counterpart; nonblank outcomes agree across all 94 shared phages | Corroborates the same-collection name mapping. Sequence/stock identity is a separate check. |
| Full public input inventory | 403 deposited host assemblies + 97 phage files; 2,066,553,341 expected bytes | A complete input audit concerns 500 files, not the 25 × 25 demo. |
| Strict file-preparation report | 104 files pass checks: eight hosts and 96 phages | The remaining 395 host files need successful retrieval. The strict parser continues to reject T4LD; its separately reviewed analysis representation is described below. No partial-set training has replaced the full experiment. |
| Published T4LD FASTA | Source bytes match the pinned Git blob; two repeated headers describe an empty record followed by 168,981 bases | Keep the original bytes and [format audit](../data/processed/genophi-source-format-audit.json). The [analysis normalization](15-relatedness-audit.md#authenticated-inputs-and-actual-searches) removes that exact-hash-reviewed empty record and preserves all actual bases. It does not rewrite the historical strict-preparer status. |
| Deposited phage sequence lengths | All 97 sources have identical first/last 60 bases; 95 of 96 metadata matches are 60 bp longer than listed, and LF110_P4 is 62 bp longer | Preserve the authenticated deposited representations. The [length diagnostic](../data/processed/genophi-phage-source-lengths.json) establishes the discrepancy; the source-generation rationale and physical termini remain unresolved. |
| Host-overlap method | skani 0.4.0 built locally; package and executable hashes recorded; policy fixed before the first host-similarity run | Ready to run once the complete host cache is available. No host-similarity output is claimed yet. |

The [local genome audit](../data/processed/genophi-genome-integrity.json) reports `complete: false`. It is a file-preparation report, not observation coverage or model accuracy. Files that pass verification are retained for resumable retrieval. The local final offline replay records missing files; the attempted online retrieval returned HTTP 403 for most host files. The requester confirmed that the tested Figshare file downloads on another computer, so the available next action is portable data retrieval.

## Reproduce preparation

Requires Python 3.11+ and curl. The metadata caches come from the [source-mapping audit](09-next-method-feasibility.md#training-source-mapping-follow-up).

```bash
uv run --project software/tipcraft python scripts/prepare-genophi-genomes.py --metadata-cache data/raw/genophi-training-map --genomes data/raw/genophi-genomes --output /tmp/genophi-genome-integrity.json --workers 4
```

Add `--offline` to authenticate the existing cache without network requests. Each download is checked against its expected size and deposited MD5 or pinned Git blob SHA-1 before atomic publication. The report adds SHA-256, contig lengths, ambiguity counts and a sequence hash preserving contig boundaries. It retains source bytes and reports repeated headers; it never merges contigs or infers new genomes. Empty/non-DNA records require review. Missing or changed files cannot become passing inputs by editing a status flag.

## Portable retrieval and compute

The official [versioned Figshare page](https://figshare.com/articles/dataset/Genome_assembly_of_the_Escherichia_Picard_collection/25941691/1) also allows downloading the host collection as an archive. Its 403 original files total about 2.06 GB before ZIP compression. Bring the downloaded ZIP back unchanged; the [archive importer](../scripts/import-genophi-hosts.py) maps original deposited filenames to the file-ID cache names automatically, verifies size/MD5/FASTA, and preserves existing verified files. It reads ZIP members without extracting their paths. Partial or unexpected contents remain explicit in its report.

```bash
uv run --project software/tipcraft python scripts/import-genophi-hosts.py --metadata-cache data/raw/genophi-training-map --archive /actual/path/to/downloaded-hosts.zip --genomes data/raw/genophi-genomes --output data/processed/genophi-host-import.json
```

Replace the archive path with the actual returned ZIP and choose a fresh report name for each import. No real archive has been imported yet; the importer has been exercised on synthetic archives and its inventory interface checked against all 403 actual metadata records. After import, replay the standard preparation audit offline, then run the frozen host-relatedness audit. Do not rename or edit the FASTA contents manually.

The [kit builder](../scripts/build-genophi-download-kit.py) packages the original downloader and authenticated public metadata. Its allowlist includes three helper scripts, their license, public source metadata, a file manifest and standalone verification instructions. It includes no private chat or target design files. Building twice produced byte-identical ZIPs, and the extracted kit's verifier and 500-file inventory were actually checked.

```bash
uv run --project software/tipcraft python scripts/build-genophi-download-kit.py --metadata-cache data/raw/genophi-training-map --output data/raw/transfer/TipCraft_genome_download_kit.zip
```

Use a new output filename if the ZIP already exists. The kit README gives standalone Python commands. Transfer the resulting `genomes/` and `genome-integrity.json` back, merge without overwriting existing verified files, then run the project preparation command with `--offline`. The known T4LD empty record can remain a reported format exception during transfer; send the other files and report rather than manually altering its sequence. The original T4LD source bytes are already retained locally for curation.

Current work needs network access and disk space, rather than a GPU. Allow about 6 GB for the source data plus transfer/archive working space. Before full protein clustering, model training or structural evaluation, measure resource use on a representative subset and record the resulting hardware requirements. The requester has offered another computer if additional compute is needed.

The [actual annotation check](16-annotation-feasibility.md) completed 110 available inputs in about 102 seconds with peak serial gene-caller RSS about 0.91 GB. This stage fits locally. It does not measure the memory required for all-host clustering or model training.

## Frozen host-relatedness triage

[The host-overlap policy](../software/tipcraft/benchmarks/retro-t7-t3/genophi-overlap-policy.json) fixes skani parameters and a conservative exclusion rule before similarity results. Exclude source hosts with estimated ANI ≥99.5% and aligned fraction ≥90% in each genome to any target reference. Also exclude named target hosts and hold borderline or missing comparisons for review. A 0.05 percentage-point ANI buffer and 0.5-point AF buffer are explicit engineering choices; they are not confidence intervals.

[skani](https://doi.org/10.1038/s41592-023-02018-3) estimates ANI and aligned fraction through sparse mapping. The [within-species genomic-unit study](https://doi.org/10.1128/mbio.02696-23) supports 99.5% as a useful conservative grouping scale; this is not an E. coli stock-identity definition. The rule is not interchangeable with the original collection's core-genome patristic threshold of 10⁻⁴. A species-level 95% ANI exclusion would remove the intended training species rather than target-near backgrounds.

Once every host file is authenticated, run from the repository root with a fresh output directory:

```bash
uv run --project software/tipcraft python scripts/audit-genophi-host-overlap.py --metadata-cache data/raw/genophi-training-map --genomes data/raw/genophi-genomes --target-genomes data/raw/wish --targets software/tipcraft/benchmarks/retro-t7-t3/wish-inputs.json --policy software/tipcraft/benchmarks/retro-t7-t3/genophi-overlap-policy.json --executable data/raw/tools/skani-0.4.0/bin/skani --output-dir data/processed/genophi-overlap/hosts-v1
```

The wrapper verifies actual inputs before running, then records the executable, policy, input and raw-matrix hashes. It has no phenotype-table input. Missing comparisons remain exclusions pending review; neither filename changes nor incomplete matrices establish independence. The [tool provenance](../data/processed/genophi-audit-tool-provenance.json) records the local skani and MMseqs2 artifacts; binaries and genome files stay in ignored raw storage.

Whole-genome phage grouping and gp17 similarity require their own pre-training specification. Relatedness exclusions, source ambiguity masks and grouped external validation must all be reflected in the final training membership. The [next-method contract](../software/tipcraft/benchmarks/retro-t7-t3/next-method-contract.json) remains the governing evaluation contract; host triage does not authorize target model inference by itself.

The [relatedness follow-up](15-relatedness-audit.md) now records actual 97-phage/two-reference searches, the exact-hash-reviewed T4LD normalization, provisional near-target exclusions and reconstructed host groups. The full target-host ANI audit still awaits the remaining host genomes.

The [source-admission follow-up](20-source-admission-and-partitions.md) verifies all 402 original-source IDs against deposit IDs, propagates exclusions through full groups and defines crossed outer blocks. Its current actual check writes HOLD only; it does not turn the eight-file host cache into an admitted training set.
