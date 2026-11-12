import numpy as np
from numba import cuda
import time


@cuda.jit
def dna_alignment_kernel(target, candidate, results):
    """Compare equal-length DNA sequences in parallel."""

    idx = cuda.grid(1)

    if idx < target.size:
        if target[idx] == candidate[idx]:
            results[idx] = 1
        else:
            results[idx] = 0


@cuda.jit
def global_alignment_kernel(
    target,
    candidate,
    matrix,
    diagonal,
    target_length,
    candidate_length,
    match_score,
    mismatch_penalty,
    gap_penalty
):
    """
    Calculate one anti-diagonal of the Needleman-Wunsch
    scoring matrix in parallel.
    """

    idx = cuda.grid(1)

    start_i = max(1, diagonal - candidate_length)
    end_i = min(target_length, diagonal - 1)

    i = start_i + idx

    if i <= end_i:

        j = diagonal - i

        if target[i - 1] == candidate[j - 1]:
            score = match_score
        else:
            score = mismatch_penalty

        diagonal_value = matrix[i - 1, j - 1] + score

        up_value = matrix[i - 1, j] + gap_penalty

        left_value = matrix[i, j - 1] + gap_penalty

        best = max(
            diagonal_value,
            up_value,
            left_value
        )

        matrix[i, j] = best


def validate_sequences(target, candidate):

    if not target or not candidate:
        raise ValueError("DNA sequences cannot be empty.")

    target = target.upper()
    candidate = candidate.upper()

    valid_bases = {"A", "T", "G", "C", "N"}

    if set(target) - valid_bases or set(candidate) - valid_bases:
        raise ValueError("Invalid DNA characters detected.")

    return target, candidate


def calculate_alignment_statistics(
    aligned_target,
    aligned_candidate,
    alignment_score,
    gpu_time
):

    matches = 0
    mismatches = []
    insertions = []
    deletions = []

    for position, (base1, base2) in enumerate(
        zip(aligned_target, aligned_candidate)
    ):

        if base1 == base2:

            if base1 != "-":
                matches += 1

        elif base1 == "-":

            insertions.append({
                "position": position,
                "target": base1,
                "candidate": base2
            })

        elif base2 == "-":

            deletions.append({
                "position": position,
                "target": base1,
                "candidate": base2
            })

        else:

            mismatches.append({
                "position": position,
                "target": base1,
                "candidate": base2
            })

    alignment_length = len(aligned_target)

    mismatch_count = (
        len(mismatches)
        + len(insertions)
        + len(deletions)
    )

    similarity = (
        matches / alignment_length * 100
        if alignment_length
        else 0
    )

    return {
        "matches": matches,
        "mismatch_count": mismatch_count,
        "mismatches": mismatches,
        "insertions": insertions,
        "deletions": deletions,
        "similarity": round(similarity, 2),
        "aligned_target": aligned_target,
        "aligned_candidate": aligned_candidate,
        "alignment_length": alignment_length,
        "alignment_score": int(alignment_score),
        "gpu_time_ms": gpu_time
    }


def compare_sequences_gpu(target, candidate):
    """
    GPU DNA comparison.

    Equal lengths:
        Uses the original fast parallel CUDA comparison.

    Different lengths:
        Uses GPU Needleman-Wunsch global alignment.
    """

    target, candidate = validate_sequences(
        target,
        candidate
    )

    target_array = np.frombuffer(
        target.encode("ascii"),
        dtype=np.uint8
    )

    candidate_array = np.frombuffer(
        candidate.encode("ascii"),
        dtype=np.uint8
    )

    # FAST PATH: Equal-length sequences

    if len(target) == len(candidate):

        d_target = cuda.to_device(target_array)

        d_candidate = cuda.to_device(candidate_array)

        d_results = cuda.device_array(
            len(target),
            dtype=np.uint8
        )

        threads_per_block = 256

        blocks_per_grid = (
            len(target) + threads_per_block - 1
        ) // threads_per_block

        start_time = time.perf_counter()

        dna_alignment_kernel[
            blocks_per_grid,
            threads_per_block
        ](
            d_target,
            d_candidate,
            d_results
        )

        results = d_results.copy_to_host()

        cuda.synchronize()

        gpu_time = (
            time.perf_counter() - start_time
        ) * 1000

        matches = int(np.sum(results))

        mismatches = []

        for position, (base1, base2) in enumerate(
            zip(target, candidate)
        ):

            if base1 != base2:

                mismatches.append({
                    "position": position,
                    "target": base1,
                    "candidate": base2
                })

        return {
            "matches": matches,
            "mismatch_count": len(mismatches),
            "mismatches": mismatches,
            "insertions": [],
            "deletions": [],
            "similarity": round(
                matches / len(target) * 100,
                2
            ),
            "aligned_target": target,
            "aligned_candidate": candidate,
            "alignment_length": len(target),
            "alignment_score": sum(
                2 if a == b else -1
                for a, b in zip(target, candidate)
            ),
            "gpu_time_ms": gpu_time
        }

    # GLOBAL ALIGNMENT PATH: Different-length sequences

    n = len(target)
    m = len(candidate)

    max_cells = 4_000_000

    if (n + 1) * (m + 1) > max_cells:
        raise ValueError(
            "GPU global alignment matrix is too large. "
            "Maximum supported matrix size is 4 million cells."
        )

    match_score = 2
    mismatch_penalty = -1
    gap_penalty = -2

    # Initialize scoring matrix

    matrix = np.zeros(
        (n + 1, m + 1),
        dtype=np.int32
    )

    for i in range(n + 1):
        matrix[i, 0] = i * gap_penalty

    for j in range(m + 1):
        matrix[0, j] = j * gap_penalty

    # Transfer data to GPU

    start_time = time.perf_counter()

    d_target = cuda.to_device(target_array)

    d_candidate = cuda.to_device(candidate_array)

    d_matrix = cuda.to_device(matrix)

    threads_per_block = 256

    # Process scoring matrix anti-diagonal by anti-diagonal

    for diagonal in range(2, n + m + 1):

        start_i = max(1, diagonal - m)

        end_i = min(n, diagonal - 1)

        cells = end_i - start_i + 1

        if cells <= 0:
            continue

        blocks_per_grid = (
            cells + threads_per_block - 1
        ) // threads_per_block

        global_alignment_kernel[
            blocks_per_grid,
            threads_per_block
        ](
            d_target,
            d_candidate,
            d_matrix,
            diagonal,
            n,
            m,
            match_score,
            mismatch_penalty,
            gap_penalty
        )

    cuda.synchronize()

    # Retrieve scoring matrix

    matrix = d_matrix.copy_to_host()

    gpu_time = (
        time.perf_counter() - start_time
    ) * 1000

    # CPU traceback

    aligned_target = []

    aligned_candidate = []

    i = n
    j = m

    while i > 0 or j > 0:

        if i > 0 and j > 0:

            score = (
                match_score
                if target[i - 1] == candidate[j - 1]
                else mismatch_penalty
            )

            diagonal_score = (
                matrix[i - 1, j - 1] + score
            )

        else:

            diagonal_score = float("-inf")

        up_score = (
            matrix[i - 1, j] + gap_penalty
            if i > 0
            else float("-inf")
        )

        left_score = (
            matrix[i, j - 1] + gap_penalty
            if j > 0
            else float("-inf")
        )

        if (
            i > 0
            and j > 0
            and matrix[i, j] == diagonal_score
        ):

            aligned_target.append(target[i - 1])

            aligned_candidate.append(candidate[j - 1])

            i -= 1
            j -= 1

        elif i > 0 and matrix[i, j] == up_score:

            aligned_target.append(target[i - 1])

            aligned_candidate.append("-")

            i -= 1

        else:

            aligned_target.append("-")

            aligned_candidate.append(candidate[j - 1])

            j -= 1

    aligned_target.reverse()

    aligned_candidate.reverse()

    aligned_target = "".join(aligned_target)

    aligned_candidate = "".join(aligned_candidate)

    return calculate_alignment_statistics(
        aligned_target,
        aligned_candidate,
        matrix[n, m],
        gpu_time
    )