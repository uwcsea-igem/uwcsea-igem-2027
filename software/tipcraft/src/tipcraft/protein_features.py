"""Train-only protein-family occurrence vocabulary and fixed projection.

This dependency-free kernel preserves the authenticated phage pilot's hit
ordering and pc_N naming. It accepts no phenotypes and runs no alignments.
Callers authenticate inputs and freeze the resulting vocabulary separately.
"""
from collections import defaultdict
from collections.abc import Mapping
import math
from numbers import Integral
import struct


RAW_SCHEMA = ["Query", "Target", "Score", "SeqIdentity", "E-value", "qStartPos",
              "qEndPos", "qLen", "tStartPos", "tEndPos", "tLen"]


def _identifier(value):
    return (isinstance(value, str) and bool(value) and value == value.strip()
            and not any(ord(c) < 32 or ord(c) == 127 for c in value))


def _ids(values, description):
    if isinstance(values, (str, bytes)):
        raise ValueError(f"Expected an iterable of {description}")
    result = list(values)
    if any(not _identifier(v) for v in result) or len(set(result)) != len(result):
        raise ValueError(f"Invalid or duplicate {description}")
    return result


def _lines(lines):
    if isinstance(lines, (str, bytes)):
        raise ValueError("Pass an iterable of text lines, not a whole text string")
    for number, line in enumerate(lines, 1):
        if not isinstance(line, str):
            raise ValueError("Expected text alignment/cluster lines")
        yield number, line.rstrip("\r\n")


def _family_mapping(families, training_ids):
    if not isinstance(families, Mapping) or set(families) != set(training_ids):
        raise ValueError("Family map must cover exactly the training protein IDs")
    for member, family in families.items():
        if not _identifier(member) or not _identifier(family) or family not in families or families[family] != family:
            raise ValueError("Unknown/invalid family representative or representative self membership")


def read_families(lines, training_protein_ids):
    """Read two-field representative/member rows, requiring complete coverage."""
    ids = set(_ids(training_protein_ids, "training protein IDs"))
    mapping = {}
    for number, line in _lines(lines):
        fields = line.split("\t")
        if len(fields) != 2:
            raise ValueError(f"Expected two cluster fields at row {number}")
        family, member = fields
        if family not in ids or member not in ids or member in mapping:
            raise ValueError("Unknown or duplicate training cluster member/representative")
        mapping[member] = family
    _family_mapping(mapping, ids)
    return mapping


def _lengths(values, description):
    if not isinstance(values, Mapping):
        raise ValueError(f"Expected {description} mapping")
    _ids(values, description)
    if any(isinstance(n, bool) or not isinstance(n, Integral) or n < 1 for n in values.values()):
        raise ValueError("Protein lengths must be positive integers")


def _float32(text):
    try:
        value = struct.unpack("f", struct.pack("f", float(text)))[0]
    except (TypeError, ValueError, OverflowError, struct.error) as exc:
        raise ValueError("Invalid float32 alignment value") from exc
    if not math.isfinite(value):
        raise ValueError("Nonfinite alignment value")
    return value


def best_hits(lines, query_lengths, training_lengths, families):
    """Stream 11-field MMseqs createtsv rows; keep raw-first float32 ties.

    Coordinates are zero-based inclusive. Targets must belong to the complete
    frozen training family map. Score does not break identity/E-value ties.
    """
    _lengths(query_lengths, "query protein lengths")
    _lengths(training_lengths, "training protein lengths")
    _family_mapping(families, training_lengths)
    best, tied, count = {}, {}, 0
    for number, line in _lines(lines):
        fields = line.split("\t")
        if len(fields) != 11:
            raise ValueError(f"Expected 11 alignment fields at row {number}")
        query, target = fields[:2]
        if query not in query_lengths or target not in training_lengths:
            raise ValueError("Alignment contains an unknown query or nontraining target")
        _, identity, evalue = map(_float32, fields[2:5])
        if not 0 <= identity <= 1 or evalue < 0:
            raise ValueError("Invalid identity/E-value")
        try:
            qs, qe, qlen, ts, te, tlen = map(int, fields[5:])
        except ValueError as exc:
            raise ValueError("Invalid alignment coordinates") from exc
        if not (qlen == query_lengths[query] and tlen == training_lengths[target]
                and 0 <= qs <= qe < qlen and 0 <= ts <= te < tlen):
            raise ValueError("Alignment coordinates/length differ from protein lengths")
        key = (-identity, evalue)
        row = {"target": target, "family": families[target], "identity_float32": identity,
               "evalue_float32": evalue, "raw_row": number}
        if query not in best or key < best[query][0]:
            best[query] = (key, row)
            tied[query] = {families[target]}
        elif key == best[query][0]:
            tied[query].add(families[target])
        count += 1
    ambiguities = [{"query": query, "families": sorted(fs), "chosen": best[query][1]}
                   for query, fs in sorted(tied.items()) if len(fs) > 1]
    return {query: row for query, (_, row) in best.items()}, {
        "raw_rows": count, "queries_with_hits": len(best),
        "queries_without_hits": len(query_lengths) - len(best),
        "cross_family_best_hit_ties": ambiguities, "raw_schema": list(RAW_SCHEMA)}


def _protein_genomes(mapping, genomes):
    if not isinstance(mapping, Mapping):
        raise ValueError("Expected protein-to-genome mapping")
    _ids(mapping, "protein IDs")
    if any(not _identifier(g) or g not in genomes for g in mapping.values()):
        raise ValueError("Protein maps to an unlisted/invalid genome")


def _hits(hits, query_ids, families):
    if not isinstance(hits, Mapping) or not set(hits) <= set(query_ids):
        raise ValueError("Hit mapping contains an unknown query")
    for hit in hits.values():
        if (not isinstance(hit, Mapping) or not _identifier(hit.get("target")) or not _identifier(hit.get("family"))
                or hit["target"] not in families
                or hit.get("family") != families[hit["target"]]):
            raise ValueError("Hit target/family is outside the frozen training family map")


def _merge_patterns(patterns):
    pattern_features, family_features = {}, {}
    for family, pattern in sorted(patterns.items()):
        key = tuple(pattern)
        if key not in pattern_features:
            pattern_features[key] = f"pc_{len(pattern_features)}"
        family_features[family] = pattern_features[key]
    return family_features


def fit_vocabulary(training_genomes, training_protein_genome, families, training_best_hits):
    """Fit singleton filtering and exact co-occurrence merging on training only.

    No-hit proteins are counted explicitly. An empty retained vocabulary is a
    valid representation result; it is not a claim that model fitting is ready.
    """
    genomes = _ids(training_genomes, "training genomes")
    if not genomes:
        raise ValueError("Training genome list must be nonempty")
    _protein_genomes(training_protein_genome, set(genomes))
    _family_mapping(families, training_protein_genome)
    _hits(training_best_hits, training_protein_genome, families)
    occurrence = defaultdict(set)
    for query, hit in training_best_hits.items():
        occurrence[hit["family"]].add(training_protein_genome[query])
    patterns = {family: [int(genome in present) for genome in genomes]
                for family, present in sorted(occurrence.items()) if len(present) > 1}
    family_features = _merge_patterns(patterns)
    features = sorted(set(family_features.values()))
    return {"schema_version": 1, "family_features": family_features,
            "retained_family_patterns": patterns, "features": features,
            "training_genomes": genomes,
            "training_protein_families": dict(sorted(families.items())),
            "training_genomes_count": len(genomes), "training_proteins": len(training_protein_genome),
            "training_proteins_with_hits": len(training_best_hits),
            "training_proteins_without_hits": len(training_protein_genome) - len(training_best_hits),
            "clustered_families": len(set(families.values())),
            "families_after_self_search": len(occurrence),
            "families_retained_more_than_one_training_genome": len(patterns),
            "features_after_exact_occurrence_consolidation": len(features)}


def _vocabulary(vocabulary):
    if not isinstance(vocabulary, Mapping) or vocabulary.get("schema_version") != 1:
        raise ValueError("Unsupported vocabulary")
    genomes = _ids(vocabulary["training_genomes"], "vocabulary training genomes")
    families = vocabulary["training_protein_families"]
    if not isinstance(families, Mapping):
        raise ValueError("Invalid frozen protein-family map")
    _family_mapping(families, families)
    patterns = vocabulary["retained_family_patterns"]
    if not isinstance(patterns, Mapping) or not set(patterns) <= set(families.values()):
        raise ValueError("Retained families outside training vocabulary")
    for pattern in patterns.values():
        if (not isinstance(pattern, list) or len(pattern) != len(genomes)
                or any(type(v) is not int or v not in (0, 1) for v in pattern) or sum(pattern) <= 1):
            raise ValueError("Invalid retained training occurrence pattern")
    expected_map = _merge_patterns(patterns)
    if (vocabulary["family_features"] != expected_map
            or vocabulary["features"] != sorted(set(expected_map.values()))):
        raise ValueError("Vocabulary feature consolidation differs from frozen training patterns")
    return families


def project(vocabulary, genomes, protein_genome, hits):
    """Project through a fixed vocabulary using binary OR, preserving zero rows."""
    families = _vocabulary(vocabulary)
    genome_ids = _ids(genomes, "projection genomes")
    _protein_genomes(protein_genome, set(genome_ids))
    _hits(hits, protein_genome, families)
    features = list(vocabulary["features"])
    found = {genome: set() for genome in genome_ids}
    for protein, hit in hits.items():
        feature = vocabulary["family_features"].get(hit["family"])
        if feature is not None:
            found[protein_genome[protein]].add(feature)
    return features, {genome: [int(feature in found[genome]) for feature in features] for genome in genome_ids}
