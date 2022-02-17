"""Tests for population module."""
import pytest
import numpy as np

import SeqSimEvo

# tests written for SeqSimEvo.Population.__init__
#                                       .__len__
#                                       .__str__
#                                       .add_change
#                                       .consensus_sequence
#                                       .copy
#                                       .delete_sequence
#                                       .get_base
#                                       .get_sample
#                                       .get_seq
#                                       .print_sample
#                                       .sample_to_string
#                                       .stats
#                                       .to_fasta
#                                       .Hamming_distance

# __init__
def test_population_init():
    """Test population init."""
    seq = SeqSimEvo.Sequence.generate_sequence(100)
    pop = SeqSimEvo.Population(seq, 10)
    assert hasattr(pop, "changed")
    assert hasattr(pop, "changes")
    assert hasattr(pop, "sequence")
    assert hasattr(pop, "n_seq")


def test_population_len():
    """Test length."""
    seq = SeqSimEvo.Sequence.generate_sequence(100)
    pop = SeqSimEvo.Population(seq, 100)
    assert len(pop) == 100
    assert pop.n_seq == 100


def test_population_str():
    """Test print string."""
    # tests that a string gets returned, not the contents
    seq = SeqSimEvo.Sequence.generate_sequence(100)
    pop = SeqSimEvo.Population(seq, 10)
    assert isinstance(
        str(pop), str
    ), f"str() doesn't return a string, but {type(str(pop))=}"


def test_population_add_change():
    """Test `add_change`."""
    seq = SeqSimEvo.Sequence.generate_sequence(100)
    pop = SeqSimEvo.Population(seq, 100)
    pop.add_change(50, 10, 0)
    assert 50 in pop.changed, "sequence 50 not in changed {}".format(pop.changed)
    assert [10, 0] in pop.get_seq(50), "{}".format(pop.changes[50])


def test_population_add_change_out_of_range():
    """Test `add_change` out of range."""
    seq = SeqSimEvo.Sequence.generate_sequence(100)
    pop = SeqSimEvo.Population(seq, 100)

    with pytest.raises(Exception) as e_info:
        pop.add_change(150, 10, 0)

    with pytest.raises(Exception) as e_info:
        pop.add_change(50, 210, 0)


def test_population_consensus_sequence():
    """Test `consensus_sequence`."""
    seq = SeqSimEvo.Sequence(seq="AAAA")
    change = np.array([[1, 2]])
    pop = SeqSimEvo.Population(seq, 20)
    for seq_id in range(10):
        pop.changes[seq_id] = change
    assert (
        str(pop.consensus_sequence()) == "ATAA"
    ), f"consensus sequence wrong. is {pop.consensus_sequence()}, should be ATAA"


def test_population_copy():
    """Test `copy`."""
    seq = SeqSimEvo.Sequence(seq="AAAAAAAAAA")
    pop = SeqSimEvo.Population(seq, 100)
    pop_copied = pop.copy()
    pop_copied.add_change(0, 0, 1)
    assert pop.changes[0] is None
    assert pop_copied.changes is not None


def test_population_get_base():
    """Test get_base."""
    seq = SeqSimEvo.Sequence(seq="AAAAAAAAAA")
    pop = SeqSimEvo.Population(seq, 100)

    assert pop.get_base(0, 1) == 0


def test_population_get_sample():
    """Test get_sample."""
    seq = SeqSimEvo.Sequence(seq="AAAAAAAAAA")
    pop = SeqSimEvo.Population(seq, 100)

    sample = pop.get_sample(10)
    assert len(sample) == 10, "sample returns wrong length"
    assert len(set(sample)) == 10, "sample returns double seqIDS"

    sample = pop.get_sample(500)
    assert len(sample) == 100, "sample returns wrong length when oversampling"
    assert len(set(sample)) == 100, "sample returns double seqIDS when oversampling"


def test_population_get_seq():
    """Test get_seq."""
    seq = SeqSimEvo.Sequence(seq="AAAAAAAAAA")
    pop = SeqSimEvo.Population(seq, 100)

    pop.add_change(0, 5, 2)

    assert [5, 2] in pop.get_seq(0)
    assert pop.get_seq(1) is None

    with pytest.raises(Exception):
        pop.get_seq(200)


def test_population_print_sample(capfd):
    """Test print_sample."""
    seq = SeqSimEvo.Sequence(seq="AAAAAAAAAA")
    pop = SeqSimEvo.Population(seq, 100)
    pop.set_changes(1, [[5, 2], [3, 1]])

    pop.print_sample([1])
    out, _ = capfd.readouterr()

    assert ("A-5-T\t1" in out) and (
        "A-3-G\t1" in out
    ), f"printed output:\n{out}\nis incomplete. should contain A-5-T\t1 and A-3-G\t1"


def test_population_get_summary():
    """Test get_summary."""
    seq = SeqSimEvo.Sequence(seq="AAAAAAAAAA")
    pop = SeqSimEvo.Population(seq, 100)
    pop.set_changes(1, [[5, 2], [3, 1]])

    string = pop.get_summary([1])
    assert ("A-5-T\t1" in string) and (
        "A-3-G\t1" in string
    ), f"returned string:\n{string}\nis incomplete. should contain A-5-T\t1 and A-3-G\t1"


def test_population_stats():
    """Test stats."""
    seq = SeqSimEvo.Sequence(seq="AAAAAAAAAA")
    pop = SeqSimEvo.Population(seq, 100)
    pop.set_changes(1, [[5, 2], [3, 1]])
    for i in range(60):
        pop.add_change(i, 6, 3)

    statistics = pop.stats()

    assert statistics["n_seq"] == 100, "N_seq statistic wrong"
    assert statistics["unmutated"] == 40, "number unmutated sequences wrong"
    assert statistics["total_mutations"] == 62, "total number of mutations is wrong"
    assert statistics["unique_mutations"] == 3, "number of unique mutations is wrong"
    assert (
        statistics["majority_mutations"] == 1
    ), "number of majority mutations is wrong"
    assert statistics["max_fraction"] == 60.0 / 100, "maximum fraction is wrong"
    assert statistics["GA_rate"] == 0, "fraction G-A is wrong"


def test_population_to_fasta():
    """Test `to_fasta`."""
    seq = SeqSimEvo.Sequence(seq="AAAAAAAAAA")
    pop = SeqSimEvo.Population(seq, 3)
    pop.set_changes(1, [[5, 2], [3, 1]])

    fasta = ">0\nAAAAAAAAAA\n>1\nAAAGATAAAA\n>2\nAAAAAAAAAA\n"
    assert pop.to_fasta([0, 1, 2]) == fasta, "fasta incorrect"


# Hamming_distance
# def test_population_hamming_distance_correct():
#     seq = SeqSimEvo.Sequence(seq="AAAAAAAAAA")
#     sim = SeqSimEvo.Simulation(n_seq_init=10, sequence=seq)
#     pop = SeqSimEvo.Population(
#         sim,
#         changes={
#             1: [[5, 2], [3, 1]],
#             2: [[5, 1], [3, 1]],
#             3: [[4, 2], [3, 1]],
#             4: [[4, 2], [3, 1]],
#             5: [[4, 2], [3, 1], [6, 1]],
#         },
#         changed=[1, 2, 3, 4, 5],
#     )
#     # action = mean
#     assert (
#         abs(pop.hamming_distance(list(range(10))) - 1.6) < 0.05
#     ), "mean hamming distance wrong"
#     # action = poisson_fit
#     assert pop.hamming_distance(list(range(10)), action="Poisson_fit") is np.nan

#     pop = SeqSimEvo.Population(sim)
#     assert pop.hamming_distance(list(range(10)), action="Poisson_fit") == 0
