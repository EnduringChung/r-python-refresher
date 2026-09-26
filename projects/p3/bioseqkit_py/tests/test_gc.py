"""Tests for bioseqkit.gc_content — done for you as the worked example."""

import pytest

from bioseqkit import gc_content


def test_percentages():
    assert gc_content("ATGC") == 50.0
    assert gc_content("GGCC") == 100.0
    assert gc_content("ATAT") == 0.0


def test_case_insensitive():
    assert gc_content("atgcatgc") == 50.0


def test_guards():
    with pytest.raises((TypeError, ValueError)):
        gc_content(["ATGC", "ATGC"])
    with pytest.raises((TypeError, ValueError)):
        gc_content(42)


def test_n_runs_and_empty():
    assert gc_content("NNNN") == 0.0
    assert gc_content("") == 0.0
