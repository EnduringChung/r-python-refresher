"""Tests for bioseqkit.kmer_counts — unskip when you implement it."""

import pytest

from bioseqkit import kmer_counts


@pytest.mark.skip(reason="unskip when you implement kmer_counts()")
def test_counts_all_kmers():
    assert kmer_counts("ATATA", 2) == {"AT": 2, "TA": 2}


@pytest.mark.skip(reason="unskip when you implement kmer_counts()")
def test_k1_is_base_composition():
    assert kmer_counts("AAGCT", 1) == {"A": 2, "C": 1, "G": 1, "T": 1}


@pytest.mark.skip(reason="unskip when you implement kmer_counts()")
def test_output_sorted_by_name():
    keys = list(kmer_counts("ATCGATCG", 3))
    assert keys == sorted(keys)


@pytest.mark.skip(reason="unskip when you implement kmer_counts()")
def test_rejects_k_too_big():
    with pytest.raises(ValueError):
        kmer_counts("AT", 3)
