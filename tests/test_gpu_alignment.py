import sys

sys.path.insert(0, "src")

import pytest

from geneweaver.gpu_alignment import compare_sequences_gpu
from geneweaver.cpu_alignment import compare_sequences


def test_gpu_availability():

    from numba import cuda

    assert cuda.is_available(), "CUDA GPU is not available."


def test_gpu_identical_sequences():

    target = "ATGCGTAC"
    candidate = "ATGCGTAC"

    result = compare_sequences_gpu(target, candidate)

    assert result["matches"] == 8
    assert result["mismatch_count"] == 0
    assert result["similarity"] == 100.0


def test_gpu_sequence_mismatch():

    target = "ATGCGTACGTAGCTAG"
    candidate = "ATGCGTTCGTAGCTAG"

    result = compare_sequences_gpu(target, candidate)

    assert result["matches"] == 15
    assert result["mismatch_count"] == 1
    assert result["similarity"] == 93.75


def test_cpu_gpu_consistency():

    target = "ATGCGTACGTAGCTAG"
    candidate = "ATGCGTTCGTAGCTAG"

    cpu_result = compare_sequences(target, candidate)

    gpu_result = compare_sequences_gpu(target, candidate)

    assert cpu_result["matches"] == gpu_result["matches"]

    assert (
        cpu_result["mismatch_count"]
        == gpu_result["mismatch_count"]
    )

    assert (
        cpu_result["alignment_score"]
        == gpu_result["alignment_score"]
    )


def test_gpu_execution_time():

    target = "ATGCGTACGTAGCTAG"
    candidate = "ATGCGTTCGTAGCTAG"

    result = compare_sequences_gpu(target, candidate)

    assert result["gpu_time_ms"] >= 0