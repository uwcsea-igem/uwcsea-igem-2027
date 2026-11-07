# Receptor dependence: source check and claim boundaries

Reviewed 2026-11-07. **The team's latest clarification distinguishes native LPS dependence from future protein-receptor redirection. That clarifies the intended questions; it does not establish receptor switching or validate the proposed control labels.** Binding, receptor dependence, entry, productive infection and direct assembly require separate evidence.

## What the team resources establish

The current Classroom stream, the complete three-item Classwork listing, the two linked laboratory-planning documents, the relevant reading-guide descriptions and accessible Drive searches were checked. The preparatory list names protein-receptor knockout categories, but supplies no completed receptor experiment or stock-level identity. The linked documents do not contain the combined table described in the new chat. The class Drive folder listing returned no children. These are bounded access/search findings, not proof that the table does not exist elsewhere. Private source URLs and retrieval notes stay under ignored `data/raw/receptor-source-audit/`.

The proposed interpretation is two questions sharing a host list: how natural T7 depends on LPS, and whether a future variant depends on a selected surface protein. This is an interpretation of the clarification, not a recovered table or a new team scope decision. No local experimental result, approved strain order or marker choice is inferred.

## Published T7 evidence changes the control assumptions

[Qimron et al., 2006, LPS Results and Table 3](https://pmc.ncbi.nlm.nih.gov/articles/PMC1748173/) directly tests T7 on Keio-background mutants. For **WT T7**, the relevant observations are:

| Host deletion | Published observation relative to WT *E. coli* | Faithful interpretation |
| --- | --- | --- |
| `waaC` | Table 3: EOP <0.1 | Upper-bounded category; not a zero-plaque label |
| `waaF` | Table 3: 0.1 ≤ EOP <0.3 | Infection remains measurable |
| `waaG` | Table 3: EOP <0.1 | Upper-bounded category; not a zero-plaque label |
| `waaL` | LPS Results: EOP similar to WT host | Not a categorical receptor-negative control |

The table's symbols encode EOP ranges. They must not be imported as binary absence/presence labels. These observations concern that publication's materials and conditions; they neither specify the team's exact stocks nor change the Ando benchmark labels. A phenotype on an altered LPS background does not identify a newly used protein receptor.

The genes also have different functions. WaaC/WaaF build inner-core heptose linkages; WaaG starts the outer-core glucose extension ([Nakao et al., 2012, Table 2 and Fig. 1](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0051241)). WaaL is the O-antigen ligase, rather than another equivalent core-building enzyme ([WaaL structure study, 2022, Abstract](https://www.nature.com/articles/s41586-022-04555-x)). The four deletions cannot be given one shared biochemical or infection label.

## What a knockout comparison establishes

Plaques on a parent and their absence on a corresponding deletion background support **conditional genetic dependence**. That comparison alone does not establish direct physical binding, receptor sufficiency, the blocked infection stage or successful fiber assembly. Keio deletions can affect neighboring gene expression, as directly measured by [Mateus et al., 2021, Results](https://journals.asm.org/doi/10.1128/msystems.00813-21). Evidence for restoration of the phenotype and independent binding/entry observations can strengthen causal interpretation; they are distinct observations to record, not assumed outcomes.

LamB, FhuA and BtuB usage by other phages motivates receptor hypotheses. It does not demonstrate that the same proteins are sufficient for T7 entry. The [existing T3/T7 display contrasts](19-compatibility-evidence-catalogue.md) already illustrate why scaffold and complementation context matter. Natural T7 LPS recognition cannot be converted into a generic protein–protein docking task.

[Baba et al., 2006, Synopsis](https://link.springer.com/article/10.1038/msb4100050) describes the original Keio mutants with a kanamycin-resistance cassette and cassette-free derivatives. Therefore marker status belongs to the exact stock record; “every Keio-derived strain is always kanR” is too broad. No cassette-removal procedure, plasmid marker or new construct is selected here.

## Consequence for TipCraft and the December gate

Preserve the ten Ando cases and WIsH failure unchanged. The new literature observations were inspected during method development: they are source-interpretation evidence, not a fresh blind test, calibration set or additional coverage. Related parent/deletion backgrounds must not later be split into nominally independent training/test cases.

Any future receptor-oriented dataset should record the exact stock and modification, receptor molecule class, phage/display context, assay endpoint, reference host, detection bounds and source figure separately. A dependence observation cannot automatically populate binding, entry or assembly fields. An interval or unmeasured endpoint stays interval-valued or unknown.

The native whole-genome GenoPHI route remains conditional on [complete source admission](20-source-admission-and-partitions.md). Its feature pilot does not encode either reference gp17 and provides no validation for knockout discrimination or new-head comparison. **HOLD for claims of protein-receptor redirection, a universal interface or assembly prediction.** The next scientific result remains an independently specified, source-authenticated prediction evaluation; these corrections prevent an apparent receptor success from being manufactured by the control definitions.
