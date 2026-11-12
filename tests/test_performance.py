import sys
import random
import time

sys.path.insert(0, "src")

from geneweaver.cpu_alignment import compare_sequences
from geneweaver.gpu_alignment import compare_sequences_gpu


def test_performance():

    print("\n===== GENEWEAVER PERFORMANCE BENCHMARK =====")

    random.seed(42)

    bases = "ATGC"
    sequence_length = 500

    target = "".join(
        random.choices(bases, k=sequence_length)
    )

    candidate = list(target)

    # Introduce exactly 10 mutations
    for _ in range(10):

        position = random.randrange(sequence_length)
        original = candidate[position]

        candidate[position] = random.choice(
            [base for base in bases if base != original]
        )

    candidate = "".join(candidate)

    print("GPU: NVIDIA RTX 2050")
    print("DNA sequence length:", sequence_length)

    # CPU benchmark
    start = time.perf_counter()

    cpu_result = compare_sequences(target, candidate)

    cpu_time = (time.perf_counter() - start) * 1000

    print("\n===== CPU RESULTS =====")
    print("Matches:", cpu_result["matches"])
    print("Mismatches:", cpu_result["mismatch_count"])
    print("CPU Time:", round(cpu_time, 4), "ms")

    # GPU warm-up
    compare_sequences_gpu(target, candidate)

    # GPU benchmark
    start = time.perf_counter()

    gpu_result = compare_sequences_gpu(target, candidate)

    gpu_total_time = (time.perf_counter() - start) * 1000

    print("\n===== GPU RESULTS =====")
    print("Matches:", gpu_result["matches"])
    print("Mismatches:", gpu_result["mismatch_count"])
    print("GPU Time:", round(gpu_total_time, 4), "ms")

    # Validate results
    assert cpu_result["matches"] == gpu_result["matches"]
    assert cpu_result["mismatch_count"] == gpu_result["mismatch_count"]

    assert (
        cpu_result["alignment_score"]
        == gpu_result["alignment_score"]
    )

    print("\n===== VALIDATION =====")
    print("SUCCESS: CPU and GPU results match!")

    speedup = cpu_time / gpu_total_time if gpu_total_time > 0 else 0

    print("\n===== PERFORMANCE COMPARISON =====")
    print("CPU Time:", round(cpu_time, 4), "ms")
    print("GPU Total Time:", round(gpu_total_time, 4), "ms")
    print("Speedup:", round(speedup, 2), "x")

    print("\nBENCHMARK COMPLETED!")