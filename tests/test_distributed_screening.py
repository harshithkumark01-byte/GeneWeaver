import sys

sys.path.insert(0, "src")

from geneweaver.distributed_alignment import distributed_screen
from geneweaver.cpu_alignment import compare_sequences


def test_parallel_candidate_screening():

    target = "ATGCGTACGTAGCTAG"

    candidates = [
        {
            "id": "Candidate_A",
            "sequence": "ATGCGTACGTAGCTAG"
        },
        {
            "id": "Candidate_B",
            "sequence": "ATGCGTTCGTAGCTAG"
        },
        {
            "id": "Candidate_C",
            "sequence": "TTTTTTTTTTTTTTTT"
        }
    ]

    result = distributed_screen(
        target,
        candidates
    )

    assert result["candidates_processed"] == 3

    assert result["execution_mode"] == (
        "Dask parallel candidate screening"
    )

    assert len(result["results"]) == 3

    # Best candidate should rank first.

    assert result["results"][0]["candidate_id"] == "Candidate_A"

    # Verify every result against CPU alignment.

    for candidate in candidates:

        expected = compare_sequences(
            target,
            candidate["sequence"]
        )

        actual = next(
            item for item in result["results"]
            if item["candidate_id"] == candidate["id"]
        )

        assert actual["alignment_score"] == expected["alignment_score"]

        assert actual["matches"] == expected["matches"]

        assert actual["similarity"] == expected["similarity"]