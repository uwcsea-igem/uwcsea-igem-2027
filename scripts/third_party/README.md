# VIRIDIC directional calculation adaptation

`viridic_directional.py` is a standalone Python adaptation of the directional
identity calculation in Cristina Moraru and colleagues' VIRIDIC. It is separate
from the MIT-licensed TipCraft package and runs through a subprocess interface.
This adaptation is licensed under GPL-3.0-only; the upstream license is included
as `VIRIDIC-LICENSE`. Python adaptation and input validation added 2026-11-07.

Primary method: Moraru et al. (2020),
[VIRIDIC—A Novel Tool to Calculate the Intergenomic Similarities of Prokaryote-Infecting Viruses](https://doi.org/10.3390/v12111268).
The reference implementation is pinned to
[`01_VIRIDIC_sim_calc.R`, commit d7db8031a27bc6b06ff5cb1245e372caf3794b93](https://github.com/CristinaMoraru/VIRIDIC/blob/d7db8031a27bc6b06ff5cb1245e372caf3794b93/VIRIDIC/stand_alone/viridic_scripts/01_VIRIDIC_sim_calc.R).

The adaptation preserves contained-HSP removal order, union-covered range
removal, directional aligned-string overlap trimming and the symmetric
whole-genome normalization. The reported original aligned-span sums retain
upstream behavior and can exceed genome length after partially overlapping
HSPs; they are not union coverage estimates. Similarities include unrounded
values and values rounded to three decimals.

Two upstream cases fail explicitly: multiple full-query HSPs enter a branch
assigning a 100000-identity sentinel, and multiple reverse-query HSPs conflict
with ascending query-coordinate ranges. Single reverse-query HSPs and reverse
subject alignments are supported. Negative or excessive identity totals also
fail rather than being clamped. Duplicate partial HSPs follow upstream's
contained-range removal order. The overlap regex is unanchored in upstream;
this behavior is preserved rather than replaced with a different algorithm.

Synthetic tests exercise these branches. The [recorded original-R comparison](../../data/processed/viridic-original-r-comparison-v2.json)
actually executes the unchanged pinned directional function on the frozen
99-genome HSP dataset: all 4,851 pair outputs and 9,702 nonself directional
counts/spans agree at the declared precision. This supports equivalence on
that dataset, not universal equivalence or reproduction of the full upstream
workflow. Historical output flags remain unchanged; consult the separate
comparison evidence. No BLAST search is launched by this utility.

```bash
python scripts/third_party/viridic_directional.py \
  --blast-tsv /absolute/path/to/blast.tsv \
  --genome-lengths /absolute/path/to/genome-lengths.json \
  --output /absolute/path/to/intergenomic.json
```

The length file is a JSON object of genome IDs and positive integer lengths.
BLAST TSV uses exactly:
`qseqid sseqid evalue bitscore qlen slen qstart qend sstart send qseq sseq nident gaps`.
The output covers every distinct unordered declared genome pair; missing
directional hits contribute zero. Preserve BLAST execution provenance separately
to distinguish a successful no-hit search from absent or failed execution.
