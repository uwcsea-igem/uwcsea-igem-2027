# Training-label provenance audit

Reviewed 2026-11-07. **Retain 37,632 source-supported binary labels and exclude 156 source-empty pairs from any future fitting, calibration or evaluation until clarified.** This is conservative curation for a future TipCraft run. No model has been trained on this subset and no target scores were produced.

## Cross-source reconciliation

The [original experiment table](https://github.com/mdmparis/coli_phage_interactions_2023/blob/99656d8c020fcd1dbbc30a02ee10274a07655562/data/interactions/interaction_matrix.csv) contains 402 hosts and 96 phage columns. The [pinned GenoPHI E. coli table](https://github.com/Noonanav/GenoPHI/blob/f7e30abc1b0bbba1ae76b8f58a3363bbbc5cd54f/data/interaction_matrices/ecoli_interaction_matrix.csv) contains 37,788 pairs across 402 hosts and 94 phages. A fixed, unique case-insensitive alphanumeric host-name correspondence reconciles punctuation changes; the phage names match exactly. Neither phenotype similarity nor target infection labels determine this mapping.

An independent reimplementation confirmed the [machine-readable label audit](../data/processed/genophi-label-audit.json):

| Source-table state for pairs present in GenoPHI | Pairs | Observed relationship |
| --- | ---: | --- |
| Nonempty source value | 37,632 | Every imported binary label equals `int(source_value > 0)`; zero disagreements |
| Empty source cell | 156 | Every imported label is 0; source cell alone does not support a negative observation |

The [audit script](../scripts/audit-genophi-labels.py) records the complete host correspondence and each source-empty pair. Its authenticated inputs and counts reproduce independently. Its `eligible_label_rows_before_sequence_exclusions` is a label-curation count, not a declaration that training-data identity, homology exclusions or evaluation independence are resolved.

## What the original raw annotations establish

The original [README, line 14](https://github.com/mdmparis/coli_phage_interactions_2023/blob/99656d8c020fcd1dbbc30a02ee10274a07655562/README.md#L14) defines raw scores as 0 for no lytic interaction, 1 for lytic interaction, and `n` for a result not interpretable in the image. It explains that these annotations precede aggregation and excludes phage-isolation hosts, which were assessed separately, from this raw file.

We inspected the [pinned raw annotation table](https://github.com/mdmparis/coli_phage_interactions_2023/blob/99656d8c020fcd1dbbc30a02ee10274a07655562/data/interactions/raw/raw_interactions.csv) for the 156 source-empty pairs:

| Raw-data finding | Count |
| --- | ---: |
| Pairs represented in this raw table | 155, each with nine annotation records |
| Represented pairs containing at least one `n` | 155 |
| Represented pairs containing at least one raw 1 as well as `n` | 18 |
| Annotation records over those 155 pairs | 1,033 `n`, 334 zeros, 28 ones |
| Source-empty pair without records in this file | One: `PDP110` / `LF73_P4` |

These records establish ambiguity, not a replacement binary outcome. In particular, the 18 pairs with a raw 1 cannot simply be relabelled positive without the experiment's aggregation criteria. Absence of the one pair from this raw file does not establish that it was untested. The complete aggregation rationale and the reason for every empty cell have not been recovered here.

The original [visualization notebook](https://github.com/mdmparis/coli_phage_interactions_2023/blob/99656d8c020fcd1dbbc30a02ee10274a07655562/data/interactions/bokeh_draw_interaction_matrix.ipynb) reads the matrix with missing values and reshapes it with `stack()`. Its [prediction script, lines 106–108](https://github.com/mdmparis/coli_phage_interactions_2023/blob/99656d8c020fcd1dbbc30a02ee10274a07655562/dev/predictions/predict_all_phages.py#L106) removes missing outcomes, but reads a local named input rather than proving the provenance of this published quantitative CSV. Neither file establishes an intended rule to interpret every source-empty cell as noninfection.

This review does **not** establish why the downstream table contains zeros for these cells, whether additional information informed that conversion, or any error by the source authors. The actionable choice is to exclude the ambiguous pairs from our own future run. Removing rows from a table also cannot retroactively remove their influence from an already trained model; a pretrained bundle would require its own training-membership audit.

## Artifact identity and replay

| Artifact | Identity |
| --- | --- |
| Original source commit | `99656d8c020fcd1dbbc30a02ee10274a07655562` |
| Original aggregate table Git blob | `25a3f894c5a3832e81b42c8de60c12ffb9429569` |
| Original aggregate table SHA-256 | `ec0aa5553a726e2161ea15d687a677ee2ed623873fb2a897432d0777582ad829` |
| GenoPHI table SHA-256 | `c4795d786a1e13ac998e263dbf4aaf95fc0f94c52b842d02b6f7fc6e99ccfef2` |
| Raw annotation Git blob, checked against pinned tree | `bee02ea5e74ad9c546bb20240c97476ac5ba6871` |
| Raw annotation SHA-256 | `56f7ba5e01d00d1780978f8922397137a42458fa0b6583c91d92425ac9bf3401` |

Aggregate inputs were read from `data/raw/genophi-training-map/original-interaction-matrix.csv` and `data/raw/genophi-audit/data/interaction_matrices/ecoli_interaction_matrix.csv`. The independent raw-data check used `/tmp/tipcraft-original-raw-interactions.csv`; raw third-party material is not committed. Reproduce from the repository root with Python 3.11+:

```bash
uv run --project software/tipcraft python scripts/audit-genophi-labels.py --source data/raw/genophi-training-map/original-interaction-matrix.csv --model-table data/raw/genophi-audit/data/interaction_matrices/ecoli_interaction_matrix.csv --output /tmp/genophi-label-audit.json
curl --fail --silent --show-error --location 'https://raw.githubusercontent.com/mdmparis/coli_phage_interactions_2023/99656d8c020fcd1dbbc30a02ee10274a07655562/data/interactions/raw/raw_interactions.csv' --output /tmp/tipcraft-original-raw-interactions.csv
uv run --project software/tipcraft python - <<'PY'
import collections
import csv
import hashlib
import json
from pathlib import Path
raw = Path('/tmp/tipcraft-original-raw-interactions.csv').read_bytes()
assert hashlib.sha256(raw).hexdigest() == '56f7ba5e01d00d1780978f8922397137a42458fa0b6583c91d92425ac9bf3401'
audit = json.loads(Path('/tmp/genophi-label-audit.json').read_text())
wanted = {(p['source_strain'], p['phage']) for p in audit['source_blank_pairs']}
by_pair = collections.defaultdict(list)
for row in csv.DictReader(raw.decode().splitlines(), delimiter=';'):
    pair = (row['bacteria'], row['phage'])
    if pair in wanted:
        by_pair[pair].append(row['score'])
assert len(wanted) == 156 and len(by_pair) == 155
assert all(len(v) == 9 and 'n' in v for v in by_pair.values())
assert wanted - by_pair.keys() == {('PDP110', 'LF73_P4')}
assert sum('1' in v for v in by_pair.values()) == 18
print(collections.Counter(s for values in by_pair.values() for s in values))
PY
```

The source references and four exchange abstentions in the T7/T3 evaluation are unchanged. The next independent scoring run still requires authenticated sequence mappings, related-group exclusions and a frozen training specification.
