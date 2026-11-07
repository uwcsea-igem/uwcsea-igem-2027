# Relatedness audit before independent prediction

Updated 2026-11-07. **Sequence inspection identifies a near-target phage hidden behind an unrelated sample name: `427_P2`.** It is present in the external training table and must be considered in the exclusion audit. This advances an independent evaluation beyond checking whether training filenames contain “T7” or “T3”. The results below describe sequence relatedness, not host-range or assembly performance.

## Authenticated inputs and actual searches

The [normalization manifest](../data/processed/genophi-phage-inputs.json) binds all 97 deposited collection sequence representations and both frozen target references to their source bytes and sequence hashes. The T7 and T3 GenBank records match the frozen whole-genome sequences; their annotated gp17 translations were independently checked against the CDS DNA. The reference proteins are `NP_042005.1` (553 aa) and `NP_523342.1` (558 aa).

The published T4LD file has a leading empty record. Normalization removes that exact-hash-reviewed empty record, retains the original source file, and preserves all 168,981 actual bases. Other real contigs are never concatenated. This resolves the format for this analysis; it does not silently change the strict downloader's historical audit.

A separate [source-length diagnostic](../data/processed/genophi-phage-source-lengths.json), governed by its [recorded policy](../software/tipcraft/benchmarks/retro-t7-t3/genophi-phage-source-length-policy.json), finds identical first/last 60 bases in **all 97** deposited sources. Of the 96 files with metadata rows, 95 are exactly 60 bp longer than the listed genome size; LF110_P4 is 62 bp longer. There are no ambiguous bases explaining the differences. `411_P3` has no metadata row. The inspected pinned README and phage directory do not establish the generation rationale. These observations do not establish physical genome termini, circularity or a justified trimming rule. The current calculations use the exact deposited representations, with no bases removed; comparisons to physical genome length require this qualification.

The [frozen phage policy](../software/tipcraft/benchmarks/retro-t7-t3/genophi-phage-policy-v2.json) specifies BLAST 2.17.0, the paper-derived nucleotide search parameters, a total-genome similarity threshold of 70%, and a separate near-gp17 rule of ≥90% amino-acid identity with ≥80% query coverage. The latter is a project exclusion policy, not a published phenotype boundary; it does not establish reciprocal full-protein coverage.

We actually ran all 99 genomes against one another and searched both reference gp17 proteins against those genomes. The completed nucleotide output contains 99 query records in two BLAST batches (67 + 32); each batch's completion count and every query's reported row count are checked. Both protein queries recover their own target reference at the declared sensitivity. The [report](../data/processed/genophi-phage-overlap.json) binds the policy, executable, normalized inputs and raw outputs; the [numeric matrix](../data/processed/genophi-phage-similarities.tsv) contains all **4,851 unordered pairs**.

## Provisional exclusions

| Collection genome | Whole-genome result involving a target | Near-gp17 result | Training-table membership |
| --- | --- | --- | --- |
| T7_Portugal | T7 total similarity 99.984% | T7 gp17: 99.458% identity, 100% query coverage | Absent from the 94-phage modelling table |
| 427_P2 | T7 total similarity 83.784%; T3 total similarity 65.453% | T3 gp17: 91.935% identity, 100% query coverage | **Present: exclude provisionally** |

The whole-genome graph places T7, T7_Portugal and 427_P2 in one target-containing component; T3 is in its own target component. The whole-genome and full-length protein checks therefore detect different aspects of relatedness. Neither result establishes an experimental donor, a phenotype, a functional module boundary or construct identity.

No other source genome meets the frozen broader gp17 diagnostic (40% identity / 80% query coverage) in this search. This does not rule out shorter shared regions or remote homologs. The full-genome and gp17 checks are grouping/exclusion tools; they are not a nearest-neighbor phenotype lookup.

## Calculation validation and source ambiguity

The whole-genome calculation is a standalone GPL-3.0-only Python adaptation of the [pinned VIRIDIC implementation](https://github.com/CristinaMoraru/VIRIDIC/blob/d7db8031a27bc6b06ff5cb1245e372caf3794b93/VIRIDIC/stand_alone/viridic_scripts/01_VIRIDIC_sim_calc.R). It remains outside the MIT TipCraft package and uses a subprocess interface. Synthetic tests cover overlap trimming and invalid input branches. Its reported span sums retain upstream behavior and must not be read as union coverage. [Adaptation and license details](../scripts/third_party/README.md).

We froze and ran two accession-specific numerical controls from the [VIRIDIC paper](https://doi.org/10.3390/v12111268). The MG969412.1/MG969413.1 calculation is **62.4455%**, consistent with the published **62.4%** at its printed precision. The accession pair NC_001604.1/NC_020483.1 gives **0.5706%**, while the paper states **0.34%** for T7/HTVC011P. Both BLAST 2.17.0 and 2.9.0 give the same calculated values in these runs; changing the version alone did not explain the difference.

Visual inspection of the actual Figure 3 confirms that its caption calls NC_020483.1 **HTVC011P** and its panel labels also use that accession. The retrieved record instead identifies it as **HTVC019P**; NCBI identifies NC_020482.1 as HTVC011P. This is an unresolved publication-name/accession correspondence. The original benchmark inputs are needed before treating the comparison as an unambiguous implementation failure or substituting another accession as a reproduced result. No name-based replacement alignment was run.

The [source-correction artifact](../data/processed/viridic-reference-source-correction.json) preserves this interpretation and binds the visually inspected primary figure and historical reports. Read it together with the [2.17 numerical control](../data/processed/viridic-reference-validation.json) and [2.9 version diagnostic](../data/processed/viridic-reference-version-validation.json). Reproduction of that publication number remains unresolved; no name-based replacement is reported as its reproduction.

### Direct comparison with the original implementation

We separately executed the **unchanged pinned original R directional function** on the already completed HSP data, using real dependency packages. The [integrity-refined comparison](../data/processed/viridic-original-r-comparison-v2.json), governed by a [pre-comparison specification](../software/tipcraft/benchmarks/retro-t7-t3/viridic-original-r-check-v2.json), passes **all 4,851 unordered pairs, 9,702 nonself directional counts/spans, and 99 self-controls**. Three-decimal outputs agree exactly; the maximum unrounded difference is **5.68 × 10⁻¹⁴ percentage points**, below the fixed 10⁻¹² tolerance. The supplied R executable/library, 27 dependency versions and 2,066 installed package files are authenticated. Missing or nonfinite numeric results fail explicitly. No BLAST search was repeated and no similarity parameter was changed.

The [first comparison](../data/processed/viridic-original-r-comparison.json) and its raw outputs remain unchanged. V2 adds executable/package checks and nonfinite-result rejection, and explicitly binds both earlier and current normalized-length files to the same hash. **The directional calculation is verified against the original function on this frozen dataset.** This does not establish universal equivalence, reproduce the entire VIRIDIC pipeline, resolve the publication's accession/name conflict or authenticate physical stock identity. The unrelated ambiguous publication number is no longer an implementation blocker for these declared sequence-grouping inputs. Final training membership still requires host exclusions and the frozen preprocessing/split specification.

## Related-host groups

The original collection's public core-genome distance matrix supports reconstructing **301 connected components across 404 source hosts** at the published strict distance threshold `<10⁻⁴`. All **402 modelling-table hosts** map into these groups. All 404 source nodes participate so a nontraining source can bridge related training members. Group IDs depend on sorted membership and remain stable under input reordering.

The [host-group report](../data/processed/genophi-host-core-groups.json) authenticates the source matrix and modelling table. BL21 is a singleton in this matrix. BW25113 and MG1655 are absent by name; that does not establish their absence by sequence similarity. These reconstructed groups are not the missing original CV fold assignments and do not replace the target-reference ANI audit. The script reads the external table for host IDs; phenotype values do not determine groups.

For the primary independent evaluation, build protein-family reference clusters and consolidated feature mappings using the training partition alone, then assign validation/target proteins to those frozen references. Learning vocabulary or centroids from held-out proteins would change the experiment into a transductive analysis; any such sensitivity study needs its own declaration and cannot replace the primary contract. Sequence use for the predeclared relatedness exclusion audit does not fit predictive features.

## Reproduction and next step

All commands run from the repository root. Prepare the frozen source caches and target genomes using the earlier audit commands. Raw published sequences and alignment strings remain in ignored storage.

```bash
uv run --project software/tipcraft python scripts/normalize-genophi-phages.py --metadata-cache data/raw/genophi-training-map --genomes data/raw/genophi-genomes --target-genomes data/raw/wish --target-manifest software/tipcraft/benchmarks/retro-t7-t3/wish-inputs.json --genbank-references data/raw/genophi-phage-audit/references --output-dir data/raw/genophi-phage-audit/normalized-v3
uv run --project software/tipcraft python scripts/audit-genophi-phage-overlap.py --policy software/tipcraft/benchmarks/retro-t7-t3/genophi-phage-policy-v2.json --normalized data/raw/genophi-phage-audit/normalized-v3 --blast-bin data/raw/tools/blast-2.17.0/ncbi-blast-2.17.0+/bin --adapter scripts/third_party/viridic_directional.py --workdir data/raw/genophi-phage-audit/alignment-v2 --output /tmp/genophi-phage-overlap.json
uv run --project software/tipcraft python scripts/group-genophi-hosts.py --matrix data/raw/genophi-training-map/original-host-core-distances.tsv --model-table data/raw/genophi-audit/data/interaction_matrices/ecoli_interaction_matrix.csv --output /tmp/genophi-host-core-groups.json
```

Normalization and a fresh alignment require new, nonexistent output directories. For the existing completed alignment directory, add `--replay` to the overlap command; this checks completion and recomputes the calculation without repeating BLAST. Raw output uses outfmt 7 so completed no-hit queries can be distinguished from missing query blocks. The recorded first alignment was stopped to add these integrity controls; the revised run retained the same alignment parameters and exclusion thresholds. A parser fix then recognized BLAST's legitimate multiple completion batches and replayed the completed search.

The original-R integrity comparison can be reproduced with the authenticated private runtime and a fresh output directory:

```bash
uv run --project software/tipcraft python scripts/compare-viridic-original-r-v2.py --policy software/tipcraft/benchmarks/retro-t7-t3/viridic-original-r-check-v2.json --r-home data/raw/tools/r-4.6.2/framework/R.framework/Versions/4.6/Resources --r-library data/raw/tools/r-4.6.2/library --runtime-manifest data/raw/tools/r-4.6.2/runtime-provenance.json --workdir data/raw/viridic-r-comparison/reproduction-1 --output /tmp/viridic-original-r-comparison.json
```

The report embeds the exact runtime/package manifest and official source URLs. It was extracted locally without a system installation; a different runtime needs its own frozen comparison specification. The helper does not install dependencies or launch alignments.

Before training: finish the host-genome transfer and target ANI exclusions; settle the [annotation correspondence or uniform-reannotation adaptation](16-annotation-feasibility.md); combine phage/host groups with the 156 source-empty label exclusions; disclose the deposited sequence representation; and freeze the complete training specification. **No new target infection scores, assembly scores or designs have been generated.** The six historical WIsH scores and four exact-construct abstentions remain unchanged.
