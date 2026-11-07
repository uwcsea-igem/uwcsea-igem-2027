"""Train-only amino-acid k-mer vocabulary and direct proteome projection.

This is an explicit global-occurrence adaptation, not the upstream GenoPHI
novel-genome protein-family assignment workflow. No phenotypes are accepted.
"""
from collections import Counter, defaultdict
from collections.abc import Mapping
from dataclasses import dataclass

AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")


def _proteomes(values):
    if not isinstance(values, Mapping) or not values:
        raise ValueError("Expected a nonempty sample-to-protein mapping")
    for sample, proteins in values.items():
        if not isinstance(sample, str) or not sample or sample != sample.strip():
            raise ValueError("Invalid sample identifier")
        if not isinstance(proteins, Mapping) or not proteins:
            raise ValueError("Every sample must supply deposited proteins")
        for protein, sequence in proteins.items():
            if not isinstance(protein, str) or not protein or protein != protein.strip():
                raise ValueError("Invalid protein identifier")
            if not isinstance(sequence, str) or not sequence or not set(sequence) <= AMINO_ACIDS:
                raise ValueError("Expected nonempty uppercase standard-amino-acid proteins")


def _kmers(sequence, k):
    return {sequence[i:i + k] for i in range(len(sequence) - k + 1)}


@dataclass(frozen=True)
class KmerVocabulary:
    """Immutable training order and merged k-mer membership.

    Members share a binary occurrence pattern in training samples. Projection
    takes their union: one member present in a new proteome activates a feature.
    This rule is explicit because members can diverge in unseen samples.
    """
    k: int
    minimum_protein_occurrences: int
    training_ids: tuple[str, ...]
    features: tuple[tuple[str, ...], ...]
    retained_kmer_count: int


def fit_vocabulary(training_proteomes, *, k=5, minimum_protein_occurrences=2):
    """Fit exclusively on training proteins; count each k-mer once per protein.

    The rare filter uses distinct proteins, including multiple proteins within
    the same sample. It does not require occurrence in multiple genomes.
    """
    for value in (k, minimum_protein_occurrences):
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            raise ValueError("k and protein-occurrence minimum must be positive integers")
    _proteomes(training_proteomes)
    ids = tuple(sorted(training_proteomes))
    protein_counts, patterns = Counter(), defaultdict(set)
    for sample in ids:
        for sequence in training_proteomes[sample].values():
            for mer in _kmers(sequence, k):
                protein_counts[mer] += 1
                patterns[mer].add(sample)
    grouped = defaultdict(list)
    retained = 0
    for mer, count in protein_counts.items():
        if count >= minimum_protein_occurrences:
            grouped[tuple(sorted(patterns[mer]))].append(mer)
            retained += 1
    features = tuple(sorted(tuple(sorted(members)) for members in grouped.values()))
    return KmerVocabulary(k, minimum_protein_occurrences, ids, features, retained)


def project_proteomes(vocabulary, proteomes):
    """Project through the frozen training vocabulary, without refitting.

    Returns sample IDs in lexical order, a matching binary matrix, and IDs of
    zero-vector samples. An all-zero vector is preserved, never dropped.
    """
    if not isinstance(vocabulary, KmerVocabulary):
        raise ValueError("Expected a fitted KmerVocabulary")
    _proteomes(proteomes)
    # The reverse index includes all equivalent training members, not just one.
    index = {mer: column for column, members in enumerate(vocabulary.features) for mer in members}
    ids, matrix, zeros = sorted(proteomes), [], []
    for sample in ids:
        present = set()
        for sequence in proteomes[sample].values():
            present.update(index[mer] for mer in _kmers(sequence, vocabulary.k) if mer in index)
        matrix.append([int(column in present) for column in range(len(vocabulary.features))])
        if not present:
            zeros.append(sample)
    return {"sample_ids": ids, "matrix": matrix, "zero_vector_sample_ids": zeros}
