# Development environment handoff

Snapshot: 2026-11-08, Asia/Singapore. Repository: https://github.com/uwcsea-igem/uwcsea-igem-2027. Branch: `codex/independent-method-feasibility`. The user paused and then cleared the persistent goal, and requested an upload and handoff; no scientific computation was resumed for this handoff.

## Current result

TipCraft is evaluating an independent strain-level prediction route on public retrospective cases. No biological interaction model has yet been fitted or scored. New-method coverage of the original ten Ando observations remains **0/10**. WIsH previously scored the six natural cases and failed the intended strain-ranking test. The four exact exchanged-gp17 construct cases still abstain because construct identity/sequences are unresolved.

The secondary BASEL–ECOR panel is frozen at 50 phages, 23 hosts and 1,093 observed eligible pairs, excluding 147 of 1,240 source observations. All 60 original unknown pairs remain unknown; 57 fall within the eligible 50-by-23 grid. Grouped 3-by-3 both-unseen folds, three model scopes and fixed untuned settings are recorded. Genotype references are not authenticated laboratory stocks. See [source-panel admission](40-secondary-panel-and-feature-launch.md), [official sequence correspondence](35-basel-phage-accession-correspondence.md), and the actual manifests rather than inferring identity from names.

V1 and v2 feature attempts exposed genuine technical control failures. V3 uses ordinary external-query projection for training and query proteins, with a training-only vocabulary frozen before external query preparation. It separately counts matching assigned families, matched absence and discrepancies. V1/v2 code and failures remain preserved. This is an adapted GenoPHI-inspired protein-family/CatBoost route, not a completed replication or validated prediction method.

## Exact execution state at handoff

The [machine-readable snapshot](../data/processed/environment-handoff-state.json) records the frozen plan hash, stale raw ledger hash, completed task IDs and migration boundary.

- **7/18 v3 feature tasks completed and were verified.** Their reports are committed under `data/processed/secondary-features-v3-*.json`. Every listed artifact hash, exact nine-command plan, normalized input and pre-query vocabulary boundary was checked.
- **2 tasks were interrupted and have no final report:** `host0-phage2-host` and `host1-phage1-host`. Their old ledger entries say running, but the producers/controller are gone and a subsequent process inspection finds no related computation still running. Partial directories and command-5 logs remain. Their overall exits are unknown, not zero.
- **9 tasks never started.** Exact IDs are in the snapshot.
- No source/native model was fitted, no partial-block accuracy was reported, and no original observation was removed.

The first accepted host partition has 816 features and a copied control of 4,960 matching assigned families, two matched absences and zero discrepancies. Its final feature names and all genome vectors equal its preserved failed v2 output; the old failure remains failed under its original contract. Bas67/Bas68 remain zero-vector phage queries in `host0-phage1` and must not be dropped. See the [seven-task checkpoint](47-v3-production-validation-checkpoint.md).

The old exec session IDs are no longer usable after the interruption. **Do not relaunch the existing batch command as a resume operation.** Its guard requires a completed initial-only validation ledger and rejects an already-started remainder. Do not overwrite that ledger, reuse partial destinations, or report an unknown process exit as zero.

## What GitHub contains and what must move separately

GitHub contains source code, tests, policies, public processed reports and documentation. Approximately **24 GB** of original inputs, tools, raw logs, intermediate databases and execution receipts live in ignored `data/raw/`; exact file-count/byte totals are in the snapshot. A clone alone cannot re-authenticate the accepted artifacts or resume the frozen execution. Copy that directory separately with an ordinary file transfer, preserve its internal structure, and verify the existing per-artifact SHA256 bindings after transfer. No raw archive, credential, private WhatsApp export or local environment is uploaded by this handoff.

Old root: `/Users/hahatt/Documents/GitHub/uwcsea-igem-2027`. Some frozen records contain absolute paths to that root; command plans also bind the original executable path. If the new root differs, preserve the original JSON bytes and hashes. Implement and document a separate audited relocation map or a new execution plan; do not bulk-rewrite old manifests and pretend they remain the same evidence.

The tool binaries are platform-specific. Do not copy virtual environments as a portable installation, assume a macOS executable runs on Linux, or replace a bound binary silently. Inspect the new OS/architecture, install the locked Python runtimes with `uv`, and version/revalidate any necessary tool-build or platform changes. The previous machine had 14 logical CPUs and about 36 GiB RAM; the v3 plan allowed two simultaneous feature tasks, four threads each. No GPU training requirement has been established.

## Verified checks

At handoff, the following commands pass: **396 package tests** and **5 locked-runtime synthetic CatBoost integration tests**. These are software/technical checks, not biological prediction accuracy. No remote CI outcome for the handoff commit is claimed.

```sh
uv run --project software/tipcraft --with biopython==1.86 python -m unittest discover -s software/tipcraft/tests
uv run --project software/tipcraft/experiments/learning --no-sync python -m unittest discover -s software/tipcraft/experiments/learning/tests
```

Install/sync the locked learning project on a fresh environment before using `--no-sync`. Check the version pins in its project/lock files.

## Next scientific actions after migration

1. Read the [copyable next-AI prompt](49-next-ai-resume-prompt.md), this snapshot, [v3 protocol](46-v3-symmetric-production-protocol.md), [latest checkpoint](47-v3-production-validation-checkpoint.md), and [gate report](07-december-gate-report.md). Verify current checkout/process/file state.
2. Transfer and authenticate the raw inputs. Resolve path/platform issues without rewriting historical evidence. Audit the two partial jobs. If recovery cannot establish completed artifacts under an explicit protocol, record the external interruption and fix a new recovery attempt with unchanged scientific settings and fresh destinations. This is not an automatic retry of a failed biological or mapping result. Preserve the seven accepted reports and all old partial artifacts.
3. Complete all 18 feature tasks under an explicit recovery plan. Retain matched absence, zero vectors, failures and exact coverage. Before the first biological fit, freeze the complete nine-block experiment manifest with actual feature hashes, exact fit commands/destinations, pinned v3 policies and full matched evaluation.
4. Execute paired, host-only and phage-only fitting/evaluation for all nine blocks. Evaluate the same eligible subset and report real scores, baselines, group weighting, informative contrast counts, failures and coverage. Never pool raw margins from different models as host-order evidence. Run the separately fixed label-free conditional-order diagnostic after authenticating scores.
5. Execute final-source/native features and the fixed native transfer fit, irrespective of disappointing source performance, provided technical prerequisites pass. Freeze all three bundles and all six native-pair scores before native label evaluation. Compare WIsH and baselines on the same six pairs and retain four construct abstentions in the ten-case denominator.
6. Perform the predeclared numeric/opaque inference replays using saved models. This checks the inference boundary, not end-to-end anonymity or blinding. Finish an English decision report with actual results and a justified GO/HOLD/NO-GO. Small retrospective samples do not establish reliable generalization.

Chat with the user in Chinese. Repository/scientific writing stays English. Use the installed anti-defensive-writing skills for clear evidence-led writing, while preserving all material negative results and limitations. Do not download the skipped Picard collection, generate new recognition heads/orderable DNA, infer assembly compatibility from host prediction, promise awards, recreate a cleared goal, or launch an automation without the user's request.

## Verified NAS migration package

The complete raw tree is now available in the NAS `AdminHome` share at `00_Inbox/TipCraft-raw-2026-11-08_225713`. On the original Mac this is `/Volumes/AdminHome/00_Inbox/TipCraft-raw-2026-11-08_225713`. The archive is `tipcraft-data-raw.tar.zst` (3,690,772,655 bytes, approximately 3.69 GB). It contains 12,567 regular files (22,376,443,674 logical bytes), 666 preserved symlinks and 1,691 directories. NAS readback verified every regular file SHA256, symlink string, directory inventory and compressed-stream integrity. Original local files remain intact.

Archive SHA256: `729139c9a71e7b3157f2de600bad9b3e4057d9f0cec7d41ea2b16e950599e8e6`. The package includes `raw-manifest.json`, `SHA256SUMS`, `verification.json`, `README-RESTORE.md`, a portable `verify-raw-tree.py`, and a copy of the next-AI prompt. Follow its staging-directory extraction instructions, then verify the restored tree before relocation or scientific execution. This archive is on the private NAS, not GitHub.
