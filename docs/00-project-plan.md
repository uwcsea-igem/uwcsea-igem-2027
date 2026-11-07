# Project plan — 2027 season

Draft for team review · 2026-11-07 · project name and gate dates unapproved.

[Long-term roadmap](10-platform-roadmap.md) · [2026 official-rule evidence for 2027 preparation](11-igem-rule-evidence.md) · [Prior art and testable innovation](12-prior-art-and-innovation.md)

**Goal.** Develop a reusable phage recognition-interface standard and TipCraft research software. The immediate deliverable is a reproducible public-case evaluation that distinguishes host prediction, binding and fiber compatibility. The long-term target-to-design workflow advances as those components acquire their own validation. The intended product is the standard, dataset and software; biological work follows verified materials, facilities, supervision and funding.

**Scope from the discussion.** On 4 November Kyungmin proposed a modular tail-fiber platform with 12–16 characterized tips, approximately 25 strains, and December/February gates. On 5 November the discussion changed to three levels: (1) natural-phage characterization, (2) standalone protein binding, (3) engineered phages. Kevin and Kyungmin favored Level 2; on 7 November Kevin retained the phage direction with a sensor backup. These are discussion records, not approvals or validated feasibility. A published swap does not establish a universal socket; standalone binding does not establish assembled-phage function.

**Deliverables before the December gate.** Four-claim literature audit with primary sources; role proposal; reproducible retrospective dataset and evaluator; documented baseline and missing predictions; assembly evidence requirements; verified facility/inventory and budget assumptions. Keep an independent sensor workspace and document all contribution evidence.

| Proposed gate | GO evidence | HOLD / NO-GO consequence |
| --- | --- | --- |
| 2026-12-15: computational feasibility | Source-reviewed labels, frozen independent scores and training exclusions, informative within-phage rankings above the declared tied/random and WIsH comparators on the same cases; report uncertainty, groups and coverage; full coverage requires all ten cases | Missing predictions, baseline-only performance or leaked labels: no claim of validated scoring; revise and retain the evidence/evaluation deliverable |
| 2026-12-15: wet-lab readiness | Named supervisor, written facility access, exact eligible material list, approved local safety procedures and funded scope including contingency | Any missing item: no wet-lab start or orders; continue computational work and assess sensor fallback |
| 2027-02-15: scope commitment | Team and supervisor review results, assembly evidence separately from binding, actual spending and remaining capacity | Narrow scope to supported outputs; choose computational and/or sensor deliverable; do not escalate automatically |

**Metrics.** Host ranking: within-phage positive–negative concordance, top-host outcome and exact observation coverage; ties and random rankings use the same candidate subset. Binary accuracy and balanced accuracy are secondary when the method supports a separately declared binary evaluation. Below-detection labels describe the source experiment. Assembly: direct-evidence coverage with supported, unsupported and unknown counts. Related literature cases provide a feasibility diagnostic; generalization needs an independent grouped holdout. Calibrate structural thresholds on external data before use.

**Budget.** Kyungmin's 5 November informal estimates were approximately S$13,000 all-new, S$7,500 realistic and S$5,000 with existing equipment/resources; these are not vendor quotes. Registration/travel and any equipment purchases need separate totals. Kevin suggested a 20% gap and approximately 15% backup allocation; treatment of those amounts remains to be agreed. See [budget](03-budget.md).

**Open decisions.** Assign exact stock IDs and endpoint-specific roles to the now-distinguished native-LPS and future protein-targeting questions; review the control assumptions in the [receptor-source check](21-receptor-dependence-and-source-check.md); confirm assembly assay definition; complete the predictor's training and exclusion specification; approve gate dates, owners, budget and project name. Suggested names only: Plug-n-Phage, TailorMade, SOCKET.

**Current evidence (7 November).** A real WIsH genome-composition run covers six native observations and abstains on four exchange cases. It ranks the documented T3 positive below both negatives. [December evidence report](07-december-gate-report.md) records raw scores, actual anonymous reruns and label checks. NO-GO for this method as the strain-ranking engine; GO for evidence/evaluation infrastructure. Assembly prediction, orderable DNA and wet-lab readiness remain unvalidated. [Literature audit](02-literature.md) · [decisions](01-decision-log.md) · [roles](05-roles.md).
