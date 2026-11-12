import sys
import random

sys.path.insert(0, "src")

from geneweaver.distributed_alignment import distributed_alignment
from geneweaver.cpu_alignment import compare_sequences


def test_distributed_alignment():

    print("\n===== GENEWEAVER DISTRIBUTED ALIGNMENT TEST =====")

    random.seed(42)

    bases = "ATGC"
    sequence_length = 500

    target = "".join(
        random.choices(bases, k=sequence_length)
    )

    candidate = list(target)

    # Introduce 10 mutations
    for _ in range(10):

        position = random.randrange(sequence_length)

        original = candidate[position]

        candidate[position] = random.choice(
            [base for base in bases if base != original]
        )

    candidate = "".join(candidate)

    # CPU alignment
    cpu_result = compare_sequences(target, candidate)

    # Dask alignment
    dask_result = distributed_alignment(
        target,
        candidate,
        chunk_size=100
    )

    print("\n===== CPU RESULTS =====")

    print("Matches:", cpu_result["matches"])
    print("Mismatches:", cpu_result["mismatch_count"])
    print("Similarity:", cpu_result["similarity"], "%")

    print("\n===== DASK RESULTS =====")

    print("Total DNA length:", dask_result["total_length"])
    print("Chunks processed:", dask_result["chunks_processed"])
    print("Matches:", dask_result["matches"])
    print("Mismatches:", dask_result["mismatch_count"])
    print("Similarity:", dask_result["similarity"], "%")

    assert cpu_result["matches"] == dask_result["matches"]
    assert cpu_result["mismatch_count"] == dask_result["mismatch_count"]
    assert cpu_result["alignment_score"] == dask_result["alignment_score"]

    print("\nSUCCESS: CPU and Dask results match!")
    print("\nDISTRIBUTED ALIGNMENT TEST COMPLETED!")