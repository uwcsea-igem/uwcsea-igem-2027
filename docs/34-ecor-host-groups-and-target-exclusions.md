# ECOR host groups and target-related exclusions

Updated 2026-11-07. **The secondary panel retains 23 candidate source-host references in 15 relatedness groups after excluding ECOR13 and ECOR25 as target-related.** All 900 directed comparisons and 30 self-controls are present; the previous five-reference control matrix is reproduced exactly. This is an input-partition result, not model admission or predictive performance.

## Frozen comparison

The [policy](../software/tipcraft/benchmarks/receptor-resource/ecor-host-relatedness-policy.json) was recorded before the first comparison of these ECOR references. It binds all [25 authenticated multi-contig ECOR inputs](33-ecor-wgs-input-authentication.md), the original three Ando-case references and the two previously authenticated author-browser references. The latter remain possible transitive bridges; no known reference is removed after seeing its relatedness results.

The executable, numerical thresholds and estimator parameters are unchanged from the [prior host policy](25-reference-host-overlap.md): skani 0.4.0, compression 125, marker compression 1,000, screen 80%, minimum reported alignment fraction 15% and one thread. A near edge uses ANI at least 99.5% and at least 90% alignment fraction in each direction, with the existing 0.05 ANI-point and 0.5 AF-point boundary buffers. Either reported direction meeting the near or buffered-boundary rule creates an undirected edge. Connected components retain all 30 reference nodes.

These are fixed project relatedness units, not taxonomic assignments, stock-identity tests or a guarantee that every relevant biological similarity has been excluded. Contigs remain separate records within each genome file; no DNA is concatenated or trimmed.

## Executed evidence

The [producer](../scripts/group-ecor-host-references.py) and [receipt](../data/processed/ecor-host-reference-groups.json) record:

| Control or result | Observed |
| --- | ---: |
| Reference inputs | 30 |
| Directed comparisons, including self | 900 / 900 |
| Failed self-controls | 0 / 30 |
| Previous five-reference directed controls | 25 / 25 exactly equal |
| Near/buffered graph edges | 27 |
| Full-reference components | 17 |
| Target-related ECOR references excluded | 2 |
| Remaining candidate ECOR references | 23 |
| Candidate ECOR components | 15 |

Missing comparisons, failed self-controls or drift in the prior 25 numeric/classification controls would leave partitioning on HOLD. None occurred. Input FASTA, dependency, executable and policy hashes were checked before and after execution; no interaction outcomes were read.

## Target-related components

The component containing the original K-12 targets comprises author Keio, original BW25113, original MG1655, ECOR13 and ECOR25. The separate original BL21/author BL21 component contains no ECOR member. Every source reference in any original-target-containing component is excluded from candidate training, including boundary and transitive links.

For the ECOR-as-reference direction:

| Source → original target | ANI (%) | AF source / target (%) | Fixed-rule classification |
| --- | ---: | ---: | --- |
| ECOR13 → BW25113 | 99.48 | 90.50 / 90.18 | Buffered boundary; exclude |
| ECOR13 → MG1655 | 99.51 | 90.73 / 90.21 | Near; exclude |
| ECOR25 → BW25113 | 99.49 | 90.18 / 92.40 | Buffered boundary; exclude |
| ECOR25 → MG1655 | 99.52 | 90.36 / 92.38 | Near; exclude |

All reverse directions and all other pairs remain in the receipt. Rounded ANI values do not establish exact sequence or stock identity. Boundary links are not silently retained because a printed value falls below 99.5.

## Grouped validation boundary

The remaining 23 ECOR references occupy 15 components. Three multi-host components must be kept intact in subsequent partitions:

- ECOR18, ECOR19, ECOR20, ECOR21;
- ECOR26, ECOR27;
- ECOR29, ECOR30, ECOR32, ECOR33, ECOR34.

The other twelve candidate references are singleton components under this rule. Group IDs are deterministic hashes of sorted full member lists. No train/test fold is selected here; there is no random-row fallback or phenotype-based movement between groups. The absence of a graph edge is not an assertion of complete evolutionary independence.

The initial [52-phage/25-host contrast counts](32-secondary-contrast-source-selection.md) remain raw-panel counts. They have not been relabelled as eligible validation counts after these exclusions. Exact phage accession reconciliation, phage-related exclusions, a complete group-aware evaluation manifest and endpoint admission still precede model fitting. Original Ando independent-method coverage remains **0/10**, with no new binding, assembly or functional-infection claim.

## Bounded phage-record route pilot

A parallel source task checked three accession-prefix representatives without selecting them by outcome. The root agent then re-parsed every saved official payload and corresponding author GenBank container using Biopython 1.86. The [pilot receipt](../data/processed/basel-phage-accession-route-pilot.json) establishes exact DNA equality across the two official formats and the existing author container for:

| Source phage | Resolved official version | Sequence length |
| --- | --- | ---: |
| Bas01 | PZ445405.1 | 50,126 |
| Bas12 | MZ501053.1 | 49,298 |
| Bas14 | NC_073062.1 | 43,627 |

PZ/MZ records were available through ENA EMBL/FASTA endpoints; the RefSeq record through NCBI EFetch GenBank/FASTA endpoints. No rotation, reverse-complement matching or sequence editing was needed. Bas01's custom record ID contains a different embedded length; the actual parsed sequence and checksum are authoritative for this comparison. The researcher had read existing source-summary prose, so this is not described as blinded research.

This pilot establishes the route for a subsequent frozen bulk reconciliation, **not authentication of all 52 phages**. The other 49 records remain unretrieved by this pilot. Exact DNA equality of public representations still does not authenticate the experimental stock.

## Replay and reproduction

The [offline replay](../scripts/replay-ecor-host-groups.py) re-authenticates dependencies and all 30 FASTAs, parses the saved raw matrix, reconstructs the complete comparison summary, checks the prior 25 controls and reproduces the full component/exclusion partition. Its [receipt](../data/processed/ecor-host-reference-groups-replay.json) confirms exact agreement without executing skani, network requests or a model.

With authenticated source caches and the pinned executable, use fresh destinations:

```bash
uv run --project software/tipcraft python scripts/group-ecor-host-references.py --workdir data/raw/genophi-external-validation/ecor-host-groups-reproduction-1 --output /tmp/ecor-host-reference-groups.json
uv run --project software/tipcraft python scripts/replay-ecor-host-groups.py --report /tmp/ecor-host-reference-groups.json --output /tmp/ecor-host-reference-groups-replay.json
```

Four new tests cover target-containing bridge components, incomplete-matrix HOLD behavior, duplicate/missing component members and deterministic group IDs. Existing tests cover complete-graph transitivity, missing directions, duplicate rows and invalid numerical comparisons. The complete 289-test suite passes with Biopython 1.86 enabled.
