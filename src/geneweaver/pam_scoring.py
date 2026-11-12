def reverse_complement(sequence):
    """
    Generate the reverse complement of a DNA sequence.
    """

    complement = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G",
        "N": "N"
    }

    return "".join(
        complement[base] for base in reversed(sequence)
    )


def find_pam_sites(sequence):
    """
    Detect SpCas9 NGG PAM sites on both DNA strands.
    """

    sequence = sequence.upper()

    valid_bases = {"A", "T", "G", "C", "N"}

    if not sequence:
        raise ValueError("DNA sequence cannot be empty.")

    if set(sequence) - valid_bases:
        raise ValueError("Invalid DNA characters detected.")

    pam_sites = []

    # Forward strand
    for position in range(len(sequence) - 2):

        pam = sequence[position:position + 3]

        if pam[1:] == "GG":

            pam_sites.append({
                "position": position,
                "pam": pam,
                "strand": "+"
            })

    # Reverse strand
    reverse_sequence = reverse_complement(sequence)

    for position in range(len(reverse_sequence) - 2):

        pam = reverse_sequence[position:position + 3]

        if pam[1:] == "GG":

            original_position = len(sequence) - position - 3

            pam_sites.append({
                "position": original_position,
                "pam": pam,
                "strand": "-"
            })

    return sorted(
        pam_sites,
        key=lambda site: (site["position"], site["strand"])
    )


def score_pam_sites(sequence):
    """
    Assign illustrative prototype scores to detected PAM sites.
    """

    pam_sites = find_pam_sites(sequence)

    for site in pam_sites:

        site["score"] = 100

        if site["pam"] == "AGG":
            site["score"] = 95

    return pam_sites