from pathlib import Path
from Bio import SeqIO


VALID_BASES = {"A", "T", "G", "C", "N"}


def parse_fasta(file_path):
    """
    Read and validate DNA sequences from a FASTA file.

    Returns:
        Dictionary containing sequence IDs and DNA sequences.
    """

    file_path = Path(file_path)

    if not file_path.exists():
        raise FileNotFoundError(
            f"FASTA file not found: {file_path}"
        )

    if not file_path.is_file():
        raise ValueError("The provided path is not a file.")

    sequences = {}

    try:
        # Automatically closes the FASTA file after reading
        with file_path.open("r", encoding="utf-8") as fasta_handle:

            records = SeqIO.parse(fasta_handle, "fasta")

            for record in records:

                sequence_id = record.id
                sequence = str(record.seq).upper()

                if not sequence:
                    raise ValueError(
                        f"Sequence '{sequence_id}' is empty."
                    )

                invalid_bases = set(sequence) - VALID_BASES

                if invalid_bases:
                    raise ValueError(
                        f"Invalid DNA bases in '{sequence_id}': "
                        f"{invalid_bases}"
                    )

                sequences[sequence_id] = sequence

    except Exception as error:
        raise ValueError(
            f"Could not parse FASTA file: {error}"
        ) from error

    if not sequences:
        raise ValueError(
            "No DNA sequences found in the FASTA file."
        )

    return sequences


def load_fasta_sequences(file_path):
    """
    Load FASTA sequences as a list of dictionaries.
    """

    sequences = parse_fasta(file_path)

    return [
        {
            "id": sequence_id,
            "sequence": sequence
        }
        for sequence_id, sequence in sequences.items()
    ]


def get_target_and_candidate(file_path):
    """
    Extract the first two sequences from a FASTA file.

    First sequence  = Target DNA
    Second sequence = Candidate DNA
    """

    sequences = parse_fasta(file_path)

    if len(sequences) < 2:
        raise ValueError(
            "The FASTA file must contain at least two DNA sequences."
        )

    sequence_items = list(sequences.items())

    target_id, target_sequence = sequence_items[0]
    candidate_id, candidate_sequence = sequence_items[1]

    if len(target_sequence) != len(candidate_sequence):
        raise ValueError(
            "Target and candidate DNA sequences must have "
            "the same length for this alignment engine. "
            f"Target length: {len(target_sequence)}, "
            f"Candidate length: {len(candidate_sequence)}"
        )

    return {
        "target_id": target_id,
        "target_sequence": target_sequence,
        "candidate_id": candidate_id,
        "candidate_sequence": candidate_sequence
    }


def chunk_sequence(sequence, chunk_size=10000):
    """
    Divide a DNA sequence into smaller chunks.
    """

    if not sequence:
        raise ValueError("DNA sequence cannot be empty.")

    if chunk_size <= 0:
        raise ValueError("Chunk size must be greater than zero.")

    return [
        sequence[index:index + chunk_size]
        for index in range(0, len(sequence), chunk_size)
    ]


def get_sequence_statistics(sequence):
    """
    Calculate basic DNA sequence statistics.
    """

    sequence = sequence.upper()

    if not sequence:
        raise ValueError("DNA sequence cannot be empty.")

    invalid_bases = set(sequence) - VALID_BASES

    if invalid_bases:
        raise ValueError(f"Invalid DNA bases: {invalid_bases}")

    length = len(sequence)

    return {
        "length": length,
        "A": sequence.count("A"),
        "T": sequence.count("T"),
        "G": sequence.count("G"),
        "C": sequence.count("C"),
        "N": sequence.count("N"),
        "GC_content": (
            (sequence.count("G") + sequence.count("C"))
            / length
        ) * 100
    }