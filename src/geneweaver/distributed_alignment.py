from dask import delayed, compute

from geneweaver.cpu_alignment import compare_sequences


def align_sequences_task(target, candidate, candidate_id):
    """
    Align one candidate DNA sequence against the target.

    Each candidate is processed as an independent Dask task.
    """

    result = compare_sequences(target, candidate)

    result["candidate_id"] = candidate_id
    result["candidate_length"] = len(candidate)

    return result


def distributed_screen(target, candidates, scheduler="threads"):
    """
    Screen multiple candidate DNA sequences using Dask.

    Each candidate alignment becomes an independent task.

    Returns:
        A dictionary containing ranked alignment results.
    """

    if not target:
        raise ValueError("Target DNA sequence cannot be empty.")

    if not candidates:
        raise ValueError("At least one candidate is required.")

    target = target.upper()

    normalized_candidates = []

    for index, candidate in enumerate(candidates):

        if isinstance(candidate, str):
            candidate_id = f"Candidate_{index + 1}"
            candidate_sequence = candidate.upper()

        elif isinstance(candidate, dict):
            candidate_id = candidate.get(
                "id",
                f"Candidate_{index + 1}"
            )

            candidate_sequence = candidate.get("sequence", "").upper()

        else:
            raise TypeError(
                "Candidates must be DNA strings or dictionaries."
            )

        if not candidate_sequence:
            raise ValueError(
                f"{candidate_id} contains an empty DNA sequence."
            )

        normalized_candidates.append(
            (candidate_id, candidate_sequence)
        )

    # Create one independent Dask task per candidate.

    tasks = []

    for candidate_id, candidate_sequence in normalized_candidates:

        task = delayed(align_sequences_task)(
            target,
            candidate_sequence,
            candidate_id
        )

        tasks.append(task)

    # Execute candidate alignments through Dask.

    results = compute(
        *tasks,
        scheduler=scheduler
    )

    # Rank candidates by similarity.

    ranked_results = sorted(
        results,
        key=lambda result: (
            result["similarity"],
            result["alignment_score"]
        ),
        reverse=True
    )

    for rank, result in enumerate(ranked_results, start=1):
        result["rank"] = rank

    return {
        "target_length": len(target),
        "candidates_processed": len(ranked_results),
        "execution_mode": "Dask parallel candidate screening",
        "scheduler": scheduler,
        "results": list(ranked_results)
    }


def distributed_alignment(target, candidate, chunk_size=10000):
    """
    Backward-compatible single-candidate alignment.

    Preserves the existing dashboard and test interface.

    Note:
        A single candidate is one complete alignment task.
        Parallelism is achieved across multiple candidates
        through distributed_screen().
    """

    if not target or not candidate:
        raise ValueError("DNA sequences cannot be empty.")

    if chunk_size <= 0:
        raise ValueError(
            "Chunk size must be greater than zero."
        )

    screening_result = distributed_screen(
        target,
        [candidate]
    )

    result = screening_result["results"][0]

    result["total_length"] = len(target)
    result["chunks_processed"] = 1
    result["execution_mode"] = "Dask global alignment"

    return result