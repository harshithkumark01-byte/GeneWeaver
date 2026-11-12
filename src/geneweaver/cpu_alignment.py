import time


def compare_sequences(target, candidate):
    """
    Perform global DNA sequence alignment using
    the Needleman-Wunsch algorithm.

    Scoring:
        Match: +2
        Mismatch: -1
        Gap: -2
    """

    if not target or not candidate:
        raise ValueError("DNA sequences cannot be empty.")

    target = target.upper()
    candidate = candidate.upper()

    valid_bases = {"A", "T", "G", "C", "N"}

    if set(target) - valid_bases or set(candidate) - valid_bases:
        raise ValueError("Invalid DNA characters detected.")

    match_score = 2
    mismatch_penalty = -1
    gap_penalty = -2

    rows = len(target) + 1
    cols = len(candidate) + 1

    # Create scoring matrix

    matrix = [
        [0] * cols
        for _ in range(rows)
    ]

    # Initialize first column

    for i in range(rows):
        matrix[i][0] = i * gap_penalty

    # Initialize first row

    for j in range(cols):
        matrix[0][j] = j * gap_penalty

    # Fill scoring matrix

    for i in range(1, rows):

        for j in range(1, cols):

            if target[i - 1] == candidate[j - 1]:

                diagonal_score = (
                    matrix[i - 1][j - 1] + match_score
                )

            else:

                diagonal_score = (
                    matrix[i - 1][j - 1] + mismatch_penalty
                )

            up_score = matrix[i - 1][j] + gap_penalty

            left_score = matrix[i][j - 1] + gap_penalty

            matrix[i][j] = max(
                diagonal_score,
                up_score,
                left_score
            )

    # Traceback

    aligned_target = []

    aligned_candidate = []

    i = len(target)
    j = len(candidate)

    while i > 0 or j > 0:

        if i > 0 and j > 0:

            if target[i - 1] == candidate[j - 1]:

                diagonal_score = (
                    matrix[i - 1][j - 1] + match_score
                )

            else:

                diagonal_score = (
                    matrix[i - 1][j - 1] + mismatch_penalty
                )

        else:

            diagonal_score = float("-inf")

        up_score = (
            matrix[i - 1][j] + gap_penalty
            if i > 0
            else float("-inf")
        )

        left_score = (
            matrix[i][j - 1] + gap_penalty
            if j > 0
            else float("-inf")
        )

        if i > 0 and j > 0 and matrix[i][j] == diagonal_score:

            aligned_target.append(target[i - 1])

            aligned_candidate.append(candidate[j - 1])

            i -= 1
            j -= 1

        elif i > 0 and matrix[i][j] == up_score:

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

    # Calculate alignment statistics

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
        if alignment_length > 0
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
        "alignment_score": matrix[-1][-1]
    }


def benchmark_cpu(target, candidate, iterations=1000):
    """
    Benchmark CPU sequence alignment performance.
    """

    if iterations <= 0:
        raise ValueError("Iterations must be greater than zero.")

    start_time = time.perf_counter()

    for _ in range(iterations):

        result = compare_sequences(
            target,
            candidate
        )

    elapsed_time = time.perf_counter() - start_time

    average_time = (
        elapsed_time / iterations
    ) * 1000

    return {
        "result": result,
        "iterations": iterations,
        "total_time_ms": elapsed_time * 1000,
        "average_time_ms": average_time
    }