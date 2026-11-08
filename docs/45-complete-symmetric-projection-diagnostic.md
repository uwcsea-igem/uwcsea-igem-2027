# Complete copied-proteome external-projection diagnostic

Updated 2026-11-08. **Ordinary external projection produces consistent outcomes for all 4,962 ECOR12 copied proteins: 4,960 matching family assignments and two matched absences.** This result supports a versioned symmetric projection procedure. It does not waive the failed v2 control or admit biological fitting.

The [prospective protocol](../software/tipcraft/benchmarks/receptor-resource/symmetric-panel-probe-policy.json), SHA256 `bfbb41ba23fa24e5b47b9406147ed50493abd0577532130e8937e9f1e1069f4a`, was fixed before the [diagnostic helper](../scripts/probe-mmseqs-symmetric-panel.py) ran. It includes the entire copied ECOR12 proteome from the failed task, rather than a favorable subset. Every source/copy amino acid string is read from the actual FASTA files and verified against the stored length and sequence hash.

The source proteins are searched from a new independent query database against the unchanged target database with the exact v2 search options. The comparison reuses the authenticated external-copy search from the failed task; that side is not independently rerun. Existing clusters, family definitions and vocabulary remain fixed. No phenotype, native label, source interaction model or native prediction is used.

## Actual results

The [completed receipt](../data/processed/mmseqs-symmetric-panel-probe.json) preserves all per-protein outcomes, fixed commands, observed process exits, raw alignment artifact bindings and existing-vocabulary vectors. All three new commands exit zero; command arguments match the frozen protocol.

| Outcome | Proteins |
| --- | ---: |
| Same assigned family in both external contexts | 4,960 |
| Unassigned in both contexts | 2 |
| One-sided absence or different assigned families | 0 |
| Total compared | 4,962 |

The projected vectors under the existing vocabulary are exactly equal. The two unassigned proteins are the same 46-aa and 81-aa proteins identified in the [v2 failure audit](44-v2-prefilter-control-failure.md). Matched absence is representation consistency, not a successful family identification. Three focused synthetic tests additionally check that matched absence cannot inflate the successful-assignment count, that one-sided absence/different families remain discrepancies, and that incomplete or duplicated comparison rosters fail.

## Consequence for the production procedure

A revised procedure should send training proteins through an independent query database, using the same ordinary projection semantics as external proteins. It must rebuild the training vocabulary from those ordinary query outcomes, retain and report unassigned proteins, and freeze the vocabulary before external query preparation. Its copied-proteome control should compare all mapping outcomes, reporting successful family matches and matched absences separately, and reject any assignment discrepancy or unequal genome vector. That is a new prospective contract, not a retroactive pass for the v2 family-identification requirement.

The current diagnostic reuses a vocabulary fitted under the older training-self procedure. It does not demonstrate the behavior of a rebuilt vocabulary, every admitted training genome, another copied proteome or the complete production pipeline. Versioned implementation/policies, complete copied-proteome production controls and the full fixed evaluation manifest must precede source/native model fitting.

Biological fitting remains on technical HOLD. All v1/v2 failures and the six accepted v2 tasks remain preserved; no partial-run accuracy is reported. New-method Ando coverage remains 0/10 and the four exact exchange cases still abstain. No binding or assembly conclusion follows.
