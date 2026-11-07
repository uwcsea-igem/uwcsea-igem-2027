# Protein-calling feasibility and correspondence

Updated 2026-11-07. **The protein-calling stage runs on the current computer, but its outputs are not an exact reproduction of GenoPHI's bundled example proteomes.** This identifies a concrete preprocessing choice to settle before fitting a model. It supplies no interaction scores and does not evaluate host specificity or assembly.

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

Choose and freeze one defensible preprocessing route before training: either establish the deposited proteomes' generation correspondence, or declare uniform reannotation of every retained training, validation and target genome as an explicit methodological adaptation. Document differences rather than mixing old example proteomes with newly called inputs. No need to optimize annotation against Ando infection outcomes.

Then complete target-host sequence exclusions, record the final curated training membership, and construct reference protein families using training partitions alone. Assign held-out and target proteins to frozen references. The [next-method contract](../software/tipcraft/benchmarks/retro-t7-t3/next-method-contract.json) still requires exact splits, features, parameters, seeds and a trained bundle before inference. Full clustering and model training need their own resource measurement; the present check establishes no full-workflow hardware requirement.
