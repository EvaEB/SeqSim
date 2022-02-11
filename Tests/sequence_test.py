"""Tests for sequence module."""
import pytest
from copy import deepcopy

import SeqSimEvo

# tests written for SeqSimEvo.Sequence.__init__
#                               .__str__
#                               .__getitem__
#                               .translate_sequence
#                               .generate_seq


# __init__
def test_seq_init_sequence_correct_length():
    sequence = SeqSimEvo.Sequence.generate_sequence(100)

    assert len(sequence.sequence) == len(sequence)


def test_seq_init_creation():
    """test for presence of the necessary variables"""
    sequence = SeqSimEvo.Sequence.generate_sequence(100)

    sequence.translation
    sequence.sequence


def test_seq_str():
    """Test string generation."""
    sequence = "ATGCTGC"
    assert str(SeqSimEvo.Sequence(seq=sequence)) == sequence


# __getitem__
def test_seq_getitem_correct():
    """Test indexing."""
    original_sequence = "ATGC"
    sequence = SeqSimEvo.Sequence(seq=original_sequence)
    assert sequence[0] == "A", f"returned {sequence[0]} when 'A' was expected"
    assert sequence[1] == "T", f"returned {sequence[1]} when 'T' was expected"
    assert sequence[2] == "G", f"returned {sequence[2]} when 'G' was expected"
    assert sequence[3] == "C", f"returned {sequence[3]} when 'C' was expected"


def test_seq_getitem_fails():
    """Test indexing out of bounds."""
    sequence = SeqSimEvo.Sequence.generate_sequence(100)
    with pytest.raises(Exception):
        _ = sequence[len(sequence) + 1]


# translate_sequence
def test_seq_translate_sequence_correct():
    """Test sequence translation to encoding."""
    original_sequence = "ATCGTTCG"
    expected_translation = [0, 2, 3, 1, 2, 2, 3, 1]
    result = SeqSimEvo.Sequence.generate_sequence(100).translate_sequence(
        original_sequence
    )
    assert list(result) == expected_translation, "{} != {}".format(
        result, expected_translation
    )


# generate_seq
def test_seq_generate_seq_correct():
    """Test generating random sequence."""
    seq_len = 123
    seq1 = SeqSimEvo.Sequence.generate_sequence(seq_len)
    assert (
        len(seq1.sequence) == seq_len
    ), f"sequence length should be 123, is {len(seq1)=}"

    seq2 = SeqSimEvo.Sequence.generate_sequence(5, base_dist=[1, 1, 1, 1])
    assert seq2.sequence == [0, 0, 0, 0, 0]


def test_seq_generate_seq_fails_neg_length():
    """Test negative input for sequence generation."""
    seq_len = -123
    with pytest.raises(Exception):
        SeqSimEvo.Sequence.generate_sequence(seq_len)


def test_seq_generate_seq_fails_non_int_length():
    """Test non integer input for sequence generation."""
    seq_len = 5.5
    with pytest.raises(Exception):
        SeqSimEvo.Sequence.generate_sequence(seq_len)
