import random
import time
import statistics
import csv
from pathlib import Path

from geneweaver.cpu_alignment import compare_sequences
from geneweaver.gpu_alignment import compare_sequences_gpu


def generate_dna(length, mutation_rate=0.001):

    random.seed(42)

    bases = ["A", "T", "G", "C"]

    target = "".join(
        random.choices(bases, k=length)
    )

    candidate = list(target)

    for i in range(length):

        if random.random() < mutation_rate:

            choices = [
                base for base in bases
                if base != target[i]
            ]

            candidate[i] = random.choice(choices)

    return target, "".join(candidate)


def measure_time(function, target, candidate, repetitions=3):

    timings = []
    result = None

    for _ in range(repetitions):

        start = time.perf_counter()

        result = function(target, candidate)

        end = time.perf_counter()

        timings.append((end - start) * 1000)

    return statistics.median(timings), result


def run_benchmark():

    sizes = [
        10_000,
        100_000,
        1_000_000,
        5_000_000
    ]

    results = []

    print("\n===== GENEWEAVER PERFORMANCE BENCHMARK =====")

    print("\nWarming up GPU...")

    compare_sequences_gpu("ATGC" * 1000, "ATGC" * 1000)

    print("GPU warm-up completed!")

    print("\n===== BENCHMARK RESULTS =====")

    for size in sizes:

        print(f"\nDNA Size: {size:,} bases")

        target, candidate = generate_dna(size)

        cpu_time, cpu_result = measure_time(
            compare_sequences,
            target,
            candidate
        )

        gpu_time, gpu_result = measure_time(
            compare_sequences_gpu,
            target,
            candidate
        )

        speedup = cpu_time / gpu_time

        validation = (
            cpu_result["matches"] == gpu_result["matches"]
            and cpu_result["mismatch_count"] == gpu_result["mismatch_count"]
        )

        print(f"CPU Time: {cpu_time:.4f} ms")
        print(f"GPU Time: {gpu_time:.4f} ms")
        print(f"GPU Speedup: {speedup:.2f}x")

        print(
            "Validation:",
            "SUCCESS" if validation else "FAILED"
        )

        results.append({
            "DNA_Size": size,
            "CPU_Time_ms": round(cpu_time, 4),
            "GPU_Time_ms": round(gpu_time, 4),
            "GPU_Speedup": round(speedup, 2),
            "Matches": gpu_result["matches"],
            "Mismatches": gpu_result["mismatch_count"],
            "Similarity_Percent": gpu_result["similarity"],
            "Validation": "SUCCESS" if validation else "FAILED"
        })

    # Save benchmark results to CSV

    output_file = Path("benchmark_results.csv")

    with output_file.open(
        "w",
        newline="",
        encoding="utf-8"
    ) as file:

        writer = csv.DictWriter(
            file,
            fieldnames=results[0].keys()
        )

        writer.writeheader()
        writer.writerows(results)

    print("\n===== CSV REPORT =====")
    print("Benchmark results saved successfully!")
    print("File:", output_file.resolve())


if __name__ == "__main__":

    run_benchmark()