import sys

sys.path.insert(0, "src")

import pytest

from geneweaver.cpu_alignment import compare_sequences


def test_identical_sequences():

    target = "ATGCGTAC"
    candidate = "ATGCGTAC"

    result = compare_sequences(target, candidate)

    assert result["matches"] == 8
    assert result["mismatch_count"] == 0
    assert result["similarity"] == 100.0


def test_sequence_mismatch():

    target = "ATGCGTACGTAGCTAG"
    candidate = "ATGCGTTCGTAGCTAG"

    result = compare_sequences(target, candidate)

    assert result["matches"] == 15
    assert result["mismatch_count"] == 1
    assert result["similarity"] == 93.75


def test_sequence_insertion():

    target = "ATGC"
    candidate = "ATGCA"

    result = compare_sequences(target, candidate)

    assert result["matches"] == 4
    assert len(result["insertions"]) == 1
    assert len(result["deletions"]) == 0

def test_alignment_score():

    target = "ATGC"
    candidate = "ATGC"

    result = compare_sequences(target, candidate)

    assert result["alignment_score"] == 8


def test_invalid_empty_sequence():

    with pytest.raises(ValueError):
        compare_sequences("", "ATGC")