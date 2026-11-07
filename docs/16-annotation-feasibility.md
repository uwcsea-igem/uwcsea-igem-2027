# Protein-calling feasibility and correspondence

Updated 2026-11-07. **All 2,980 deposited example proteins translate exactly from their declared intervals in the authenticated deposited DNA.** Our fixed gene-calling procedure selects different coordinate sets and is now declared as a uniform-reannotation adaptation for the next experiment. The original 0/25 exact-proteome result remains unchanged. No interaction scores have been produced.

## Fixed check and actual inputs

The [check specification](../software/tipcraft/benchmarks/retro-t7-t3/genophi-annotation-check.json) was recorded before annotation or example comparisons. It fixes Prodigal 2.6.3 for hosts (`-p single -g 11`) and Prodigal-gv 2.11.0 for phages (`-p meta`), serial execution, an unchanged input inventory and an exact protein-multiset comparison. These versions occur in the [GenoPHI paper's software list](https://doi.org/10.1038/s41564-026-02482-5); the [Prodigal-gv implementation](https://github.com/apcamargo/prodigal-gv/tree/80e488b3cdd7123684f001759eaf50cdb879c81a) documents meta-mode use. This is a fixed candidate annotation procedure, not evidence that it reproduces the authors' exact annotation commands.

Actual scope: all **97 deposited phage representations**, both frozen phage references, all **eight locally verified external host assemblies**, and the three frozen host references: **110 genomes**. The eight available hosts are an availability subset used to check this stage; they do not substitute for the 403-file host inventory or define a training set. Source bytes, contig boundaries and order are retained; input headers become neutral identifiers. Host single-mode fits its coding model to each genome independently and uses no infection labels. No cross-genome predictive features are learned.

The [source-length diagnostic](15-relatedness-audit.md#authenticated-inputs-and-actual-searches) remains relevant: deposited phage files share terminal 60-base overlap and differ from metadata lengths. This run annotates the authenticated deposited representations without trimming them or inferring physical termini. Exact source hashes and gene-caller commands are in the [result manifest](../data/processed/genophi-annotation-check.json).

## Observed outputs

| Check | Result | Meaning |
| --- | --- | --- |
| Annotation completion | 110/110 inputs; 61,239 protein calls | This frozen procedure executes on the inspected inputs. Counts are predicted genes, not experimentally verified proteins. |
| Native target gp17 recovery | T7 NP_042005.1 and T3 NP_523342.1 each recovered once with exact amino-acid sequence hashes | The reference gp17 sequences survive this procedure. This is not host-range or binding validation. |
| Pinned example correspondence | 0/25 complete protein multisets match; 2,842 exact protein copies shared with 2,980 deposited copies and 2,975 new calls | The procedure is not an exact reproduction of these example proteomes. No alternative tool or parameter was selected to improve this comparison. |
| Serial resource observation | 102.03 seconds; peak child RSS 911,163,392 bytes (about 0.91 GB) | Adequate local capacity for this annotation check. The RSS covers serial gene-caller child processes, not subsequent clustering or model training. |

The exact comparison ignores record names, order, wrapping and one terminal stop symbol, while retaining duplicate protein copies. It rejects internal stops and empty records. All 25 deposited example files are authenticated against their pinned Git blobs; they are technical controls, with no phenotype-table input. The shared-copy count is a proteome-correspondence statistic, **not prediction accuracy**. No observed mismatch is automatically attributed to a tool bug, version difference, sequence error or biological difference.

## DNA correspondence and source follow-up

A [separately frozen CDS check](../software/tipcraft/benchmarks/retro-t7-t3/genophi-demo-cds-check.json) uses the deposited headers' one-based inclusive coordinates, strand, partial flags and code 11. It translates those exact intervals without searching other frames, shifting coordinates, trimming DNA or rerunning a gene caller. Initiation and partial-record handling follow the [Prodigal output documentation](https://github.com/hyattpd/Prodigal/wiki/Understanding-the-Prodigal-Output) and [NCBI code 11](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG11). The [result](../data/processed/genophi-demo-cds-audit.json) passes **2,980/2,980** exact translation/start/stop checks, including the seven partial records. Independent review also checked their actual contig-edge positions; the script itself checks declared partial flags without enforcing physical-edge proximity.

The deposited and regenerated sets share **2,842 exact coordinate triples**. On each side, **79** other calls share a strand and stop position with the other set but select a different start. A further **59 deposited** and **54 regenerated** calls have no same-strand stop counterpart. These are descriptions of coordinate differences, not an attribution to a particular software version or input-processing cause. Sixteen calls in each set extend beyond the metadata-listed genome length, so the length discrepancy alone does not justify removing the final 60 bases. The check preserves the original evidence and establishes sequence correspondence for the annotated intervals; it does not validate biological gene boundaries or experimental-stock identity.

The [source-provenance record](../data/processed/genophi-annotation-source-provenance.json) authenticates the pinned notebook. Zero-based **notebook cell 18** (code-cell ordinal 13), source lines 11 and 21, selects existing host/phage protein directories; lines 32–43 copy `.faa` files into the demo. Cell 10 uses the same input directories. The phage directory name includes `prodigal-gv`; it does not establish a version, mode, flags or exact DNA input. The inspected code consumes pre-called proteins and does not supply their exact generation command. Notebook prediction outputs were not inspected for this follow-up.

Reproduce the CDS audit after the earlier annotation check:

```bash
uv run --project software/tipcraft python scripts/audit-genophi-demo-cds.py --policy software/tipcraft/benchmarks/retro-t7-t3/genophi-demo-cds-check.json --output /tmp/genophi-demo-cds-audit.json
```

Use a new output path. The check authenticates the deposited examples, DNA, previous annotation outputs and declared policies before comparison. No source sequences or gene calls are changed.

## Reproduction

Run from the repository root after preparing the frozen sequence caches and all 25 pinned `data/test_data/phage_AAs/*.faa` files. They belong to [GenoPHI commit f7e30abc](https://github.com/Noonanav/GenoPHI/tree/f7e30abc1b0bbba1ae76b8f58a3363bbbc5cd54f), whose tree the check authenticates. Original `.faa`, `.gff`, logs and normalized input files remain under ignored `data/raw/genophi-annotation/run-1/`.

The local tool builds used these exact source commits and commands, without system installation:

```bash
git clone --filter=blob:none --no-checkout https://github.com/hyattpd/Prodigal.git data/raw/tools/prodigal-2.6.3/source
git -C data/raw/tools/prodigal-2.6.3/source checkout --detach 004218fbeb2cf54fe12d88ec2eede6a27d70ebf7
make -C data/raw/tools/prodigal-2.6.3/source -j 2
git clone --filter=blob:none --no-checkout https://github.com/apcamargo/prodigal-gv.git data/raw/tools/prodigal-gv-2.11.0/source
git -C data/raw/tools/prodigal-gv-2.11.0/source checkout --detach 80e488b3cdd7123684f001759eaf50cdb879c81a
make -C data/raw/tools/prodigal-gv-2.11.0/source -j 2
uv run --project software/tipcraft python scripts/check-genophi-annotation.py --policy software/tipcraft/benchmarks/retro-t7-t3/genophi-annotation-check.json --workdir data/raw/genophi-annotation/reproduction-1 --output /tmp/genophi-annotation-check.json
```

The check requires a new output directory and verifies the declared executable, provenance and source hashes before launching annotation. The existing specification records this ARM64 build, including compiler and executable hashes in its bound [Prodigal](../data/processed/genophi-prodigal-build.json) and [Prodigal-gv](../data/processed/genophi-prodigal-gv-build.json) provenance records. Those versioned JSON copies preserve the exact bytes of the corresponding raw `provenance.json` files. A rebuild on another platform requires a separately recorded tool-provenance record and check specification before execution; do not replace the frozen hashes or present a different binary as the original run. Prodigal tools retain their GPL licenses and stay outside the MIT package.

## Consequence for the next independent model

The [preprocessing component policy](../software/tipcraft/benchmarks/retro-t7-t3/genophi-preprocessing-policy.json) now fixes **uniform reannotation** with the same pinned callers and commands for all retained training, validation and target genomes. Deposited example proteins remain technical controls and are not mixed with new calls as model inputs. Keep contig boundaries and the declared deposited sequence representations; later protein-feature input removes one terminal stop per protein while retaining sequence, length, multiplicity and counted ambiguous residues. These choices were fixed before feature clustering, model fitting and new target scores; they were not selected using Ando infection outcomes. This is a named adaptation, not exact replication of upstream preprocessing.

The gene callers' distributed coding models are used as supplied; their training-genome overlap has not been audited here. They are coding models rather than supervised infection predictors, so do not describe this as generic pretraining independence. Host single mode fits a coding model to each genome independently. Neither operation reads infection labels.

Then complete target-host sequence exclusions, record the final curated training membership, and construct reference protein families using training partitions alone. Assign held-out and target proteins to frozen references. The [next-method contract](../software/tipcraft/benchmarks/retro-t7-t3/next-method-contract.json) still requires exact splits, features, parameters, seeds and a trained bundle before inference. Full clustering and model training need their own resource measurement; the present check establishes no full-workflow hardware requirement.
