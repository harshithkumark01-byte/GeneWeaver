import sys

sys.path.insert(0, "src")

from geneweaver.parser import (
    load_fasta_sequences,
    chunk_sequence,
    get_sequence_statistics
)


def test_fasta_parser():

    records = load_fasta_sequences("data/sample_genome.fasta")

    print("\n===== GENEWEAVER FASTA PARSER TEST =====")

    print("Total sequences:", len(records))

    for record in records:

        sequence = record["sequence"]

        stats = get_sequence_statistics(sequence)

        print(
            f"{record['id']} - "
            f"Length: {stats['length']} bases"
        )

    print("\nFirst sequence chunks:")

    chunks = chunk_sequence(records[0]["sequence"], 10)

    for chunk in chunks:
        print(chunk)

    print("\nPARSER TEST COMPLETED SUCCESSFULLY!")

    assert len(records) > 0
    assert all("id" in record for record in records)
    assert all("sequence" in record for record in records)


def test_chunk_sequence():

    sequence = "ATGCGTACGT"

    chunks = chunk_sequence(sequence, 5)

    assert chunks == ["ATGCG", "TACGT"]


def test_sequence_statistics():

    stats = get_sequence_statistics("ATGC")

    assert stats["length"] == 4
    assert stats["A"] == 1
    assert stats["T"] == 1
    assert stats["G"] == 1
    assert stats["C"] == 1
    assert stats["GC_content"] == 50.0