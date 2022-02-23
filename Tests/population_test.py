"""Tests for population module."""
import pytest
import numpy as np

from SeqSimEvo import Population, Sequence


def test_population_init():
    """Test population init."""
    seq = Sequence.generate_sequence(100)
    pop = Population(seq, 10)
    assert hasattr(pop, "n_seq"), "Missing attribute `n_seq`."
    assert hasattr(pop, "haplotypes"), "Missing attribute `haplotypes`."
    assert hasattr(pop, "sequence"), "Missing attribute `sequence`."
    assert hasattr(pop, "changed"), "Missing attribute `changed`."


def test_population_len():
    """Test length."""
    seq = Sequence.generate_sequence(100)
    pop = Population(seq, 100)
    assert len(pop) == 100
    assert pop.n_seq == 100


def test_population_str():
    """Test print string."""
    # tests that a string gets returned, not the contents
    seq = Sequence.generate_sequence(100)
    pop = Population(seq, 10)
    assert isinstance(
        str(pop), str
    ), f"str() doesn't return a string, but {type(str(pop))=}"


def test_population_merge():
    """Test `Population.merge`."""
    seq1 = Sequence.generate_sequence(100)
    seq2 = Sequence.generate_sequence(100)
    pop1 = Population(seq1, 10)
    pop2 = Population(seq1, 10)
    pop3 = Population(seq2, 10)
    merged = Population.merge(pop1, pop2)
    assert len(merged) == 20, f"Merged population with incorrect length {len(merged)}"
    assert merged.sequence == seq1, "Merged population with incorrect sequence."
    with pytest.raises(ValueError):
        _ = Population.merge(pop1, pop2, pop3)


def test_population_add():
    """Test `Population.__add__`."""
    seq1 = Sequence.generate_sequence(100)
    seq2 = Sequence.generate_sequence(100)
    pop1 = Population(seq1, 10)
    pop2 = Population(seq1, 10)
    pop3 = Population(seq2, 10)
    added = pop1 + pop2
    assert len(added) == 20, f"Merged population with incorrect length {len(added)}"
    assert added.sequence == seq1, "Merged population with incorrect sequence."
    with pytest.raises(ValueError):
        _ = pop1 + pop3


def test_population_get_haplotype():
    """Test `Population.get_haplotypes`."""
    seq = Sequence.generate_sequence(100)
    pop = Population(seq, 10)
    for seq_id in range(10):
        haplotype = pop.get_haplotype(seq_id)
        assert haplotype.sequence == seq, f"Incorrect haplotype sequence for {seq_id=}."
    with pytest.raises(IndexError):
        _ = pop.get_haplotype(-1)
    with pytest.raises(IndexError):
        _ = pop.get_haplotype(10)


def test_population_get_base():
    """Test `Population.get_base`."""
    seq = Sequence.generate_sequence(100)
    pop = Population(seq, 10)
    base = pop.get_base(0, 0)
    assert base == seq.sequence[0], f"Found {base} instead of {seq.sequence[0]}."


def test_population_get_sample():
    """Test get_sample."""
    seq = Sequence(seq="AAAAAAAAAA")
    pop = Population(seq, 100)

    sample = pop.get_sample(10)
    assert len(sample) == 10, "sample returns wrong length"
    assert len(set(sample)) == 10, "sample returns double seqIDS"

    sample = pop.get_sample(500)
    assert len(sample) == 100, "sample returns wrong length when oversampling"
    assert len(set(sample)) == 100, "sample returns double seqIDS when oversampling"


def test_population_add_change():
    """Test `add_change`."""
    seq = Sequence.generate_sequence(100)
    pop = Population(seq, 100)
    pop.add_change(50, 10, 0)
    assert 50 in pop.changed, f"sequence 50 not in changed {pop.changed}"
    assert [10, 0] in pop.get_haplotype(50).changes, f"{pop.get_haplotype(50).changes}"


def test_population_add_change_out_of_range():
    """Test `add_change` out of range."""
    seq = Sequence.generate_sequence(100)
    pop = Population(seq, 100)
    with pytest.raises(Exception):
        pop.add_change(150, 10, 0)
    with pytest.raises(Exception):
        pop.add_change(50, 210, 0)


def test_population_consensus_sequence():
    """Test `consensus_sequence`."""
    seq = Sequence(seq="AAAA")
    pop = Population(seq, 20)
    for seq_id in range(10):
        pop.add_change(seq_id, 1, 2)
    assert (
        str(pop.consensus_sequence()) == "ATAA"
    ), f"consensus sequence wrong. is {pop.consensus_sequence()}, should be ATAA"


def test_population_copy():
    """Test `copy`."""
    seq = Sequence(seq="AAAAAAAAAA")
    pop = Population(seq, 100)
    pop_copied = pop.copy()
    pop_copied.add_change(0, 0, 1)
    assert pop.haplotypes[0].changes is None
    assert pop_copied.haplotypes[0].changes is not None


def test_population_print_sample(capfd):
    """Test print_sample."""
    seq = Sequence(seq="AAAAAAAAAA")
    pop = Population(seq, 100)
    pop.add_change(1, 5, 2)
    pop.add_change(1, 3, 1)

    pop.print_sample([1])
    out, _ = capfd.readouterr()

    assert ("A-5-T\t1" in out) and (
        "A-3-G\t1" in out
    ), f"printed output:\n{out}\nis incomplete. should contain A-5-T\t1 and A-3-G\t1"


def test_population_get_summary():
    """Test get_summary."""
    seq = Sequence(seq="AAAAAAAAAA")
    pop = Population(seq, 100)
    pop.add_change(1, 5, 2)
    pop.add_change(1, 3, 1)

    string = pop.get_summary([1])
    assert ("A-5-T\t1" in string) and (
        "A-3-G\t1" in string
    ), f"returned string:\n{string}\nis incomplete. should contain A-5-T\t1 and A-3-G\t1"


def test_population_stats():
    """Test stats."""
    seq = Sequence(seq="AAAAAAAAAA")
    pop = Population(seq, 100)
    pop.add_change(1, 5, 2)
    pop.add_change(1, 3, 1)

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
    seq = Sequence(seq="AAAAAAAAAA")
    pop = Population(seq, 3)
    pop.add_change(1, 5, 2)
    pop.add_change(1, 3, 1)

    fasta = ">0\nAAAAAAAAAA\n>1\nAAAGATAAAA\n>2\nAAAAAAAAAA\n"
    assert pop.to_fasta([0, 1, 2]) == fasta, "fasta incorrect"


# Hamming_distance
# def test_population_hamming_distance_correct():
#     seq = Sequence(seq="AAAAAAAAAA")
#     sim = SeqSimEvo.Simulation(n_seq_init=10, sequence=seq)
#     pop = Population(
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

#     pop = Population(sim)
#     assert pop.hamming_distance(list(range(10)), action="Poisson_fit") == 0
