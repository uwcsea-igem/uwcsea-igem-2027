# SPDX-License-Identifier: GPL-3.0-only
# Python adaptation of VIRIDIC by Cristina Moraru et al.; modified 2026-11-07.
# Upstream source and GPL notice: see README.md and VIRIDIC-LICENSE alongside.
"""Adapt VIRIDIC's directional identity calculation to validated BLAST TSV.

This is a Python adaptation, not verified upstream-equivalent. It preserves
the pinned R function's HSP ordering and aligned-string trimming, including
its original-span diagnostic. Unsupported upstream sentinel/error branches
raise exceptions rather than supply apparently valid similarities.
"""
import argparse
from dataclasses import dataclass
import json
import math
import os
from pathlib import Path
import re
import tempfile
from typing import Iterable, Mapping


UPSTREAM_COMMIT = "d7db8031a27bc6b06ff5cb1245e372caf3794b93"
UPSTREAM_SOURCE = (
    "https://github.com/CristinaMoraru/VIRIDIC/blob/" + UPSTREAM_COMMIT
    + "/VIRIDIC/stand_alone/viridic_scripts/01_VIRIDIC_sim_calc.R"
)
METHOD_LABEL = "Python adaptation of pinned VIRIDIC directional identity; upstream equivalence not verified"


class UnsupportedUpstreamCase(ValueError):
    """The R function enters an unsupported or nonsensical edge-case branch."""


def _integer(value, field):
    if type(value) is not int:
        raise ValueError(f"{field} must be an integer")


@dataclass(frozen=True)
class HSP:
    qseqid: str
    sseqid: str
    evalue: float
    bitscore: float
    qlen: int
    slen: int
    qstart: int
    qend: int
    sstart: int
    send: int
    qseq: str
    sseq: str
    nident: int
    gaps: int

    def __post_init__(self):
        for name in ("qseqid", "sseqid"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value or any(x.isspace() for x in value):
                raise ValueError(f"Invalid {name}")
        for name in ("evalue", "bitscore"):
            value = getattr(self, name)
            if not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
                raise ValueError(f"{name} must be finite and nonnegative")
        for name in ("qlen", "slen", "qstart", "qend", "sstart", "send", "nident", "gaps"):
            _integer(getattr(self, name), name)
        for prefix in ("q", "s"):
            length = getattr(self, prefix + "len")
            start, end = getattr(self, prefix + "start"), getattr(self, prefix + "end")
            if length <= 0 or not 1 <= start <= length or not 1 <= end <= length:
                raise ValueError("BLAST coordinates outside genome length")
            sequence = getattr(self, prefix + "seq")
            if not isinstance(sequence, str):
                raise ValueError("BLAST aligned sequences must be strings")
            sequence = sequence.upper()
            if not sequence or not set(sequence) <= set("ACGTRYSWKMBDHVN-"):
                raise ValueError("Invalid BLAST aligned DNA string")
            if len(sequence.replace("-", "")) != abs(end - start) + 1:
                raise ValueError("BLAST coordinates disagree with aligned sequence span")
            object.__setattr__(self, prefix + "seq", sequence)
        if len(self.qseq) != len(self.sseq) or any(q == s == "-" for q, s in zip(self.qseq, self.sseq)):
            raise ValueError("Invalid BLAST alignment columns")
        if self.gaps != self.qseq.count("-") + self.sseq.count("-"):
            raise ValueError("BLAST gaps disagree with aligned strings")
        if self.nident != sum(q == s and q != "-" for q, s in zip(self.qseq, self.sseq)):
            raise ValueError("BLAST nident disagrees with aligned strings")


def _validate_lengths(lengths: Mapping[str, int]) -> None:
    for identifier, length in lengths.items():
        if not isinstance(identifier, str) or not identifier or any(x.isspace() for x in identifier):
            raise ValueError("Invalid genome identifier")
        _integer(length, "genome length")
        if length <= 0:
            raise ValueError("Genome lengths must be positive")


def parse_blast_tsv(text: str, genome_lengths: Mapping[str, int] | None = None) -> list[HSP]:
    """Read upstream's exact 14-column outfmt; retain row order and duplicates."""
    if genome_lengths is not None:
        _validate_lengths(genome_lengths)
    observed_lengths = {}
    result = []
    for number, line in enumerate(text.splitlines(), 1):
        if not line.strip():
            continue
        columns = line.split("\t")
        if len(columns) != 14:
            raise ValueError(f"BLAST row {number}: expected 14 tab-separated columns")
        try:
            hsp = HSP(columns[0], columns[1], float(columns[2]), float(columns[3]),
                      *(int(value) for value in columns[4:10]), columns[10], columns[11],
                      int(columns[12]), int(columns[13]))
            for identifier, length in ((hsp.qseqid, hsp.qlen), (hsp.sseqid, hsp.slen)):
                if identifier in observed_lengths and observed_lengths[identifier] != length:
                    raise ValueError("Inconsistent BLAST genome lengths")
                observed_lengths[identifier] = length
                if genome_lengths is not None and genome_lengths.get(identifier) != length:
                    raise ValueError("BLAST genome ID/length differs from declared inputs")
        except (ValueError, OverflowError) as exc:
            raise ValueError(f"BLAST row {number}: {exc}") from exc
        result.append(hsp)
    return result


@dataclass(frozen=True)
class DirectionalCounts:
    nident: int
    original_aligned_query_bases: int


def _overlaps(rows: list[HSP]) -> bool:
    last_end = 0
    for row in sorted(rows, key=lambda hsp: hsp.qstart):
        if row.qstart <= last_end:
            return True
        last_end = max(last_end, row.qend)
    return False


def _covered_by_other_rows(row: HSP, others: list[HSP]) -> bool:
    """Equivalent to min(IRanges coverage across this row) > 1."""
    cursor = row.qstart
    for other in sorted(others, key=lambda hsp: hsp.qstart):
        if other.qend < cursor:
            continue
        if other.qstart > cursor:
            return False
        cursor = max(cursor, other.qend + 1)
        if cursor > row.qend:
            return True
    return False


def directional_nident(hsps: Iterable[HSP]) -> DirectionalCounts:
    """Adapt pinned nident_fun, without treating its span sum as union coverage.

    Single HSPs retain upstream's abs(query-end - query-start) span. Multiple
    reverse-query HSPs cannot enter upstream's ascending IRanges machinery.
    Duplicate full-query HSPs enter upstream's 100000 sentinel branch and are
    rejected here; this is intentionally not a silent replacement algorithm.
    """
    rows = list(hsps)
    if not rows:
        return DirectionalCounts(0, 0)
    first = rows[0]
    if any((h.qseqid, h.sseqid, h.qlen, h.slen) !=
           (first.qseqid, first.sseqid, first.qlen, first.slen) for h in rows):
        raise ValueError("Directional HSPs must belong to one genome pair")
    if len(rows) == 1:
        return DirectionalCounts(first.nident, abs(first.qend - first.qstart) + 1)
    full_query = [h for h in rows if h.qstart == 1 and h.qend == first.qlen]
    if len(full_query) == 1:
        return DirectionalCounts(full_query[0].nident, first.qlen)
    if len(full_query) > 1:
        raise UnsupportedUpstreamCase("Multiple full-query HSPs enter upstream's 100000 identity sentinel branch")
    if any(h.qstart > h.qend for h in rows):
        raise UnsupportedUpstreamCase("Multiple reverse-query HSPs are unsupported by upstream IRanges ordering")
    if not _overlaps(rows):
        return DirectionalCounts(sum(h.nident for h in rows), sum(h.qend - h.qstart + 1 for h in rows))

    # Stable descending query-span order, matching dplyr::arrange ties.
    rows.sort(key=lambda h: h.qend - h.qstart + 1, reverse=True)
    for a in range(len(rows) - 1, 0, -1):
        for b in range(a):
            if rows[b].qstart <= rows[a].qstart and rows[a].qend <= rows[b].qend:
                del rows[a]
                break
    if _overlaps(rows):
        # Upstream removes ranges covered entirely by the union of others,
        # in reverse current-row order, recomputing coverage after each removal.
        for i in range(len(rows) - 1, -1, -1):
            if _covered_by_other_rows(rows[i], rows[:i] + rows[i + 1:]):
                del rows[i]
    identities = [h.nident for h in rows]
    if _overlaps(rows):
        rows.sort(key=lambda h: h.qstart, reverse=True)
        identities = [h.nident for h in rows]
        ends = [h.qend for h in rows]
        for a in range(len(rows) - 1, 0, -1):
            for b in range(a - 1, -1, -1):
                if ends[a] >= rows[b].qstart and rows[a].qstart < rows[b].qstart:
                    overlap = ends[a] - rows[b].qstart + 1
                    suffix = rows[a].qseq.replace("-", "")[-overlap:]
                    # R str_extract is deliberately unanchored and selects
                    # the first match, even if the query suffix is repeated.
                    pattern = "".join(re.escape(base) + "-*" for base in suffix)
                    match = re.search(pattern, rows[a].qseq)
                    if match is None:
                        raise UnsupportedUpstreamCase("Upstream overlap regex produced no aligned-query match")
                    query_part = match.group()
                    subject_part = rows[a].sseq[-len(query_part):]
                    removed = sum(q == s for q, s in zip(query_part, subject_part))
                    identities[a] = rows[a].nident - removed
                    ends[a] = rows[b].qstart - 1
    count = sum(identities)
    if any(value < 0 for value in identities) or count > first.qlen:
        raise UnsupportedUpstreamCase("Upstream trimming produces invalid directional identity counts")
    return DirectionalCounts(count, sum(h.qend - h.qstart + 1 for h in rows))


def calculate_similarities(hsps: Iterable[HSP], genome_lengths: Mapping[str, int]) -> list[dict]:
    """Return all distinct unordered pairs, including declared no-hit pairs."""
    _validate_lengths(genome_lengths)
    grouped = {}
    for hsp in hsps:
        if genome_lengths.get(hsp.qseqid) != hsp.qlen or genome_lengths.get(hsp.sseqid) != hsp.slen:
            raise ValueError("BLAST genome ID/length differs from declared inputs")
        grouped.setdefault((hsp.qseqid, hsp.sseqid), []).append(hsp)
    directions = {pair: directional_nident(rows) for pair, rows in grouped.items()}
    ids = sorted(genome_lengths)
    result = []
    for index, query in enumerate(ids):
        for subject in ids[index + 1:]:
            forward = directions.get((query, subject), DirectionalCounts(0, 0))
            reverse = directions.get((subject, query), DirectionalCounts(0, 0))
            similarity = 100 * (forward.nident + reverse.nident) / (genome_lengths[query] + genome_lengths[subject])
            result.append({
                "query": query, "subject": subject,
                "query_length": genome_lengths[query], "subject_length": genome_lengths[subject],
                "query_nident": forward.nident, "subject_nident": reverse.nident,
                "query_original_aligned_bases": forward.original_aligned_query_bases,
                "subject_original_aligned_bases": reverse.original_aligned_query_bases,
                "similarity_percent_unrounded": similarity, "similarity_percent": round(similarity, 3),
            })
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--blast-tsv", type=Path, required=True)
    parser.add_argument("--genome-lengths", type=Path, required=True,
                        help="JSON object mapping genome IDs to positive integer lengths")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        for source in (args.blast_tsv, args.genome_lengths):
            if args.output.resolve() == source.resolve() or (args.output.exists() and args.output.samefile(source)):
                raise ValueError("Output must be separate from inputs")
        lengths = json.loads(args.genome_lengths.read_text())
        if not isinstance(lengths, dict):
            raise ValueError("Genome lengths must be a JSON object")
        pairs = calculate_similarities(parse_blast_tsv(args.blast_tsv.read_text(), lengths), lengths)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    report = {"schema_version": 1, "method": METHOD_LABEL, "upstream_commit": UPSTREAM_COMMIT,
              "upstream_source": UPSTREAM_SOURCE, "license": "GPL-3.0-only",
              "upstream_equivalence_verified": False, "pairs": pairs}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    handle, temporary = tempfile.mkstemp(prefix=args.output.name + ".", suffix=".part", dir=args.output.parent)
    try:
        with os.fdopen(handle, "w") as stream:
            stream.write(json.dumps(report, indent=2, sort_keys=True) + "\n")
        Path(temporary).replace(args.output)
    finally:
        Path(temporary).unlink(missing_ok=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
