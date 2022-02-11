"""Tests for population module."""
import pytest
from copy import deepcopy

import SeqSimEvo

# tests written for SeqSimEvo.Population.__init__
#                                       .__len__
#                                       .__str__
#                                       .add_change
#                                       .add_sequence
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
def test_population_init_CORRECT():
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


def test_population_add_sequence():
    """Test `add_sequence`."""
    seq = SeqSimEvo.Sequence.generate_sequence(100)
    pop = SeqSimEvo.Population(seq, 100)

    new_change = [[10, 1], [20, 2]]
    new = pop.add_sequence(changes=new_change)
    assert new == 100, f"new sequence has wrong ID assigned. is {new}, should be 100"
    assert ([10, 1] in pop.changes[new]) & (
        [20, 2] in pop.changes[new]
    ), f"changes weren't added correctly. is {pop.changes[new]}, should be {new_change}"


def test_population_consensus_sequence():
    """Test `consensus_sequence`."""
    seq = SeqSimEvo.Sequence(seq="AAAA")
    pop = SeqSimEvo.Population(seq, 10)
    for _ in range(10):
        pop.add_sequence(changes=[[1, 2]])
    assert (
        str(pop.consensus_sequence()) == "ATAA"
    ), "consensus sequence wrong. is {}, should be {}".format(
        pop.consensus_sequence(), "ATAA"
    )


def test_population_copy_CORRECT():
    """Test `copy`."""
    seq = SeqSimEvo.Sequence(seq="AAAAAAAAAA")
    pop = SeqSimEvo.Population(seq, 100)

    for i in range(10):
        pop.add_sequence(changes=[[i, 2]])

    pop_copied = pop.copy()

    for i in range(10):
        pop_copied.add_sequence(changes=[[i, 3]])

    assert len(pop_copied) != len(pop)
    assert pop_copied.changes != pop.changes


# delete_sequence
def test_population_delete_sequence_CORRECT():
    seq = SeqSimEvo.Sequence(seq="AAAAAAAAAA")
    pop = SeqSimEvo.Population(seq, 1)

    id1 = pop.add_sequence(changes=[[2, 2]])
    id2 = pop.add_sequence(changes=[[3, 3]])

    pop.delete_sequence(id1)

    assert [3, 3] in pop.get_seq(id1), "{}".format(id2)
    assert pop.n_seq == 2


# get_base
def test_population_get_base_CORRECT():
    seq = SeqSimEvo.Sequence(seq="AAAAAAAAAA")
    pop = SeqSimEvo.Population(seq, 100)

    assert pop.get_base(0, 1) == 0


# get_sample
def test_population_get_sample_CORRECT():
    seq = SeqSimEvo.Sequence(seq="AAAAAAAAAA")
    pop = SeqSimEvo.Population(seq, 100)

    sample = pop.get_sample(10)
    assert len(sample) == 10, "sample returns wrong length"
    assert len(set(sample)) == 10, "sample returns double seqIDS"

    sample = pop.get_sample(500)
    assert len(sample) == 100, "sample returns wrong length when oversampling"
    assert len(set(sample)) == 100, "sample returns double seqIDS when oversampling"


# get_seq
def test_population_get_seq():
    seq = SeqSimEvo.Sequence(seq="AAAAAAAAAA")
    pop = SeqSimEvo.Population(seq, 100)

    _id = pop.add_sequence(changes=[[5, 2]])

    assert [5, 2] in pop.get_seq(_id)
    assert pop.get_seq(0) is None

    with pytest.raises(Exception):
        pop.get_seq(200)


# print_sample
def test_population_print_sample(capfd):
    seq = SeqSimEvo.Sequence(seq="AAAAAAAAAA")
    pop = SeqSimEvo.Population(seq, 100, changes={1: [[5, 2], [3, 1]]}, changed=[1])

    pop.print_sample([1])
    out, _ = capfd.readouterr()

    assert ("A-5-T\t1" in out) and (
        "A-3-G\t1" in out
    ), "printed output:\n{}\nis incomplete. should contain A-5-T\t1 and A-3-G\t1".format(
        out
    )

    with pytest.raises(Exception):
        pop.print_sample([200])


# sample_to_string
def test_population_sample_to_string():
    seq = SeqSimEvo.Sequence(seq="AAAAAAAAAA")
    pop = SeqSimEvo.Population(seq, 100, changes={1: [[5, 2], [3, 1]]}, changed=[1])

    string = pop.sample_to_string([1])
    assert ("A-5-T\t1" in string) and (
        "A-3-G\t1" in string
    ), "returned string:\n{}\nis incomplete. should contain A-5-T\t1 and A-3-G\t1".format(
        string
    )


# stats
def test_population_stats():
    seq = SeqSimEvo.Sequence(seq="AAAAAAAAAA")
    pop = SeqSimEvo.Population(seq, 100, changes={1: [[5, 2], [3, 1]]}, changed=[1])
    for i in range(60):
        pop.add_change(i, 6, 3)
    pop.add_sequence(changes=[[5, 2], [1, 1]])

    statistics = pop.stats()

    assert statistics["n_seq"] == 101, "N_seq statistic wrong"
    assert statistics["unmutated"] == 40, "number unmutated sequences wrong"
    assert statistics["total_mutations"] == 64, "total number of mutations is wrong"
    assert statistics["unique_mutations"] == 4, "number of unique mutations is wrong"
    assert (
        statistics["majority_mutations"] == 1
    ), "number of majority mutations is wrong"
    assert statistics["max_fraction"] == 60.0 / 101, "maximum fraction is wrong"
    assert statistics["GA_rate"] == 0, "fraction G-A is wrong"


def test_population_to_fasta():
    """Test `to_fasta`."""
    seq = SeqSimEvo.Sequence(seq="AAAAAAAAAA")
    pop = SeqSimEvo.Population(seq, 3, changes={1: [[5, 2], [3, 1]]}, changed=[1])

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
