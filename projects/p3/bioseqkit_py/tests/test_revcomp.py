"""Tests for bioseqkit.rev_comp — unskip when you implement it."""

import pytest

from bioseqkit import rev_comp


@pytest.mark.skip(reason="unskip when you implement rev_comp()")
def test_reverses_and_complements():
    assert rev_comp("ATGC") == "GCAT"
    assert rev_comp("A") == "T"
    assert rev_comp("AAAAAC") == "GTTTTT"


@pytest.mark.skip(reason="unskip when you implement rev_comp()")
def test_case_insensitive():
    assert rev_comp("atgc") == "GCAT"


@pytest.mark.skip(reason="unskip when you implement rev_comp()")
def test_handles_n():
    assert rev_comp("ATGN") == "NCAT"


@pytest.mark.skip(reason="unskip when you implement rev_comp()")
def test_guards():
    with pytest.raises((TypeError, ValueError)):
        rev_comp(["AT", "GC"])
