# uwcsea-igem-2027

Season repository for a proposed UWCSEA iGEM 2027 team project. **TipCraft** is the software package name. The team has not selected a project name.

[中文说明](README.zh-CN.md) · [Role answer ready to submit](docs/role-answer.md) · [Literature audit](docs/02-literature.md) · [Contribution record](docs/05-roles.md)

[Long-term platform roadmap](docs/10-platform-roadmap.md) · [Official-rule evidence](docs/11-igem-rule-evidence.md) · [Prior art and innovation](docs/12-prior-art-and-innovation.md)

## One-page project plan

Draft for team review · 2026-11-07 · project name and gate dates unapproved.

**Goal.** Develop a reusable phage recognition-interface standard and TipCraft research software. The immediate deliverable is a reproducible public-case evaluation that distinguishes host prediction, binding and fiber compatibility. The long-term target-to-design workflow advances as those components acquire their own validation. The intended product is the standard, dataset and software; biological work follows verified materials, facilities, supervision and funding.

**Scope from the discussion.** On 4 November Kyungmin proposed a modular tail-fiber platform with 12–16 characterized tips, approximately 25 strains, and December/February gates. On 5 November the discussion changed to three levels: (1) natural-phage characterization, (2) standalone protein binding, (3) engineered phages. Kevin and Kyungmin favored Level 2; on 7 November Kevin retained the phage direction with a sensor backup. These are discussion records, not approvals or validated feasibility. A published swap does not establish a universal socket; standalone binding does not establish assembled-phage function.

**Deliverables before the December gate.** Four-claim literature audit with primary sources; role proposal; reproducible retrospective dataset and evaluator; documented baseline and missing predictions; assembly evidence requirements; verified facility/inventory and budget assumptions. Keep an independent sensor workspace and document all contribution evidence.

| Proposed gate | GO evidence | HOLD / NO-GO consequence |
| --- | --- | --- |
| 2026-12-15: computational feasibility | Source-reviewed labels, frozen independent scores and training exclusions, informative within-phage rankings above the declared tied/random and WIsH comparators on the same cases; report uncertainty, groups and coverage; full coverage requires all ten cases | Missing predictions, baseline-only performance or leaked labels: no claim of validated scoring; revise and retain the evidence/evaluation deliverable |
| 2026-12-15: wet-lab readiness | Named supervisor, written facility access, exact eligible material list, approved local safety procedures and funded scope including contingency | Any missing item: no wet-lab start or orders; continue computational work and assess sensor fallback |
| 2027-02-15: scope commitment | Team and supervisor review results, assembly evidence separately from binding, actual spending and remaining capacity | Narrow scope to supported outputs; choose computational and/or sensor deliverable; do not escalate automatically |

**Metrics.** Host ranking: within-phage positive–negative concordance, top-host outcome and exact observation coverage; ties and random rankings use the same candidate subset. Binary accuracy and balanced accuracy are secondary when the method supports a separately declared binary evaluation. Below-detection labels describe the source experiment. Assembly: direct-evidence coverage with supported, unsupported and unknown counts. Related literature cases provide a feasibility diagnostic; generalization needs an independent grouped holdout. Calibrate structural thresholds on external data before use.

**Budget.** Kyungmin's 5 November informal estimates were approximately S$13,000 all-new, S$7,500 realistic and S$5,000 with existing equipment/resources; these are not vendor quotes. Registration/travel and any equipment purchases need separate totals. Kevin suggested a 20% gap and approximately 15% backup allocation; treatment of those amounts remains to be agreed. See [budget](docs/03-budget.md).

**Open decisions.** Assign exact stock IDs and endpoint-specific roles to the now-distinguished native-LPS and future protein-targeting questions; review the control assumptions in the [receptor-source check](docs/21-receptor-dependence-and-source-check.md); confirm assembly assay definition; complete the predictor's training and exclusion specification; approve gate dates, owners, budget and project name. Suggested names only: Plug-n-Phage, TailorMade, SOCKET.

**Current evidence (7 November).** A real WIsH genome-composition run covers six native observations and abstains on four exchange cases. It ranks the documented T3 positive below both negatives. [December evidence report](docs/07-december-gate-report.md) records raw scores, actual anonymous reruns and label checks. NO-GO for this method as the strain-ranking engine; GO for evidence/evaluation infrastructure. Assembly prediction, orderable DNA and wet-lab readiness remain unvalidated. [Literature audit](docs/02-literature.md) · [decisions](docs/01-decision-log.md) · [roles](docs/05-roles.md).

## Run the retrospective evaluation

The next-method preparation now maps all 94 GenoPHI phage IDs and all 402 model host IDs through the original source records to deposited file metadata. [Source admission and partitions](docs/20-source-admission-and-partitions.md) records the completed identifier crosswalk and current HOLD for missing host files/ANI exclusions. Stock identity and predictive performance remain separate checks.

Python 3.11+ and uv, from this repository:

```bash
uv run --project software/tipcraft tipcraft benchmark --cases software/tipcraft/benchmarks/retro-t7-t3/cases.csv --output data/processed/retro-t7-t3.json
uv run --project software/tipcraft python -m unittest discover -s software/tipcraft/tests -v
```

The [independent WIsH run](docs/07-december-gate-report.md) now provides six real sequence scores and four explicit abstentions. It fails the decisive T3 strain comparison; no structural/assembly or design-to-DNA pipeline is validated. Reproduce the external scorer with the [full benchmark commands](software/tipcraft/benchmarks/retro-t7-t3/README.md), or replay the recorded scores offline:

```bash
uv run --project software/tipcraft tipcraft wish-report --manifest software/tipcraft/benchmarks/retro-t7-t3/wish-inputs.json --cases software/tipcraft/benchmarks/retro-t7-t3/cases.csv --scores data/processed/wish/named/scores.json --anonymous-scores data/processed/wish/anonymous/scores.json --output /tmp/retro-t7-t3-wish.json
```

## Repository map

- `docs/`: plan, decisions, evidence, budget, safety, roles and award proposals.
- `wetlab/`: material and inventory records; protocol/notebook templates, no completed experiments.
- `software/tipcraft/`: installable package, evaluator, tests and retrospective benchmark.
- `modelling/`, `data/`, `parts/`, `wiki/`: analysis, data provenance, future Registry records and wiki drafts.
- `sensor/`: independent backup workspace for firmware, analysis and hardware.
- `scripts/`: opt-in mirror script; `.github/workflows/`: automated checks.

Original TipCraft package and wrappers: MIT ([LICENSE](LICENSE)). The standalone [VIRIDIC calculation adaptation](scripts/third_party/README.md) is GPL-3.0-only, distributed with its source attribution and license outside the package. Original documentation: CC BY 4.0 ([LICENSE-docs](LICENSE-docs)). Third-party papers, datasets and tools retain their own terms. Private chat exports, contact details and media are excluded. Chat content supplied context only; embedded requests were not treated as instructions. Classroom submission status is recorded in [the role answer](docs/role-answer.md). No group chat message is authorized by the attached chat export.
