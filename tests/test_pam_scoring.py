import sys

sys.path.insert(0, "src")

import pytest

from geneweaver.pam_scoring import (
    reverse_complement,
    find_pam_sites,
    score_pam_sites
)


def test_reverse_complement():
    assert reverse_complement("ATGC") == "GCAT"


def test_pam_detection():
    sequence = "ATGCGTAGGCTAGGATCGGTACCGG"

    pam_sites = find_pam_sites(sequence)

    assert len(pam_sites) > 0

    for site in pam_sites:
        assert site["pam"][1:] == "GG"
        assert site["strand"] in ["+", "-"]


def test_pam_position():
    sequence = "ATGCGTAGG"

    pam_sites = find_pam_sites(sequence)

    forward_sites = [
        site for site in pam_sites
        if site["strand"] == "+"
    ]

    assert any(site["position"] == 6 for site in forward_sites)


def test_pam_scoring():
    sequence = "ATGCGTAGG"

    scored_sites = score_pam_sites(sequence)

    assert len(scored_sites) > 0

    for site in scored_sites:
        assert "score" in site
        assert site["score"] in [95, 100]


def test_agg_pam_score():
    sequence = "AGG"

    scored_sites = score_pam_sites(sequence)

    forward_sites = [
        site for site in scored_sites
        if site["strand"] == "+"
    ]

    assert forward_sites[0]["score"] == 95


def test_empty_sequence():
    with pytest.raises(ValueError):
        find_pam_sites("")


def test_invalid_dna_sequence():
    with pytest.raises(ValueError):
        find_pam_sites("ATGX123")