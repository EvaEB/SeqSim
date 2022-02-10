"""Simulation module tests."""

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
def test_Population_init_CORRECT():
    seq = SeqSimEvo.Sequence.generate_sequence(100)
    pop = SeqSimEvo.Population(seq, 10)
    assert hasattr(pop, "changed")
    assert hasattr(pop, "changes")
    assert hasattr(pop, "sequence")
    assert hasattr(pop, "n_seq")


# __len__
def test_Population_len_CORRECT():
    seq = SeqSimEvo.Sequence.generate_sequence(100)
    pop = SeqSimEvo.Population(seq, 100)
    assert len(pop) == 100
    assert pop.n_seq == 100


# __str__
def test_Population_str_CORRECT():
    # tests that a string gets returned, not the contents
    seq = SeqSimEvo.Sequence.generate_sequence(100)
    pop = SeqSimEvo.Population(seq, 10)
    assert isinstance(
        str(pop), str
    ), f"str() doesn't return a string, but {type(str(pop))=}"


# add_change
def test_Population_add_change_CORRECT():
    seq = SeqSimEvo.Sequence.generate_sequence(100)
    pop = SeqSimEvo.Population(seq, 100)
    pop.add_change(50, 10, 0)
    assert 50 in pop.changed, "sequence 50 not in changed {}".format(pop.changed)
    assert [10, 0] in pop.get_seq(50), "{}".format(pop.changes[50])


def test_Population_add_change_FAILS_outofrange():
    seq = SeqSimEvo.Sequence.generate_sequence(100)
    pop = SeqSimEvo.Population(seq, 100)

    with pytest.raises(Exception) as e_info:
        pop.add_change(150, 10, 0)

    with pytest.raises(Exception) as e_info:
        pop.add_change(50, 210, 0)


# add_sequence
def test_Population_add_sequence_CORRECT():
    seq = SeqSimEvo.Sequence.generate_sequence(100)
    pop = SeqSimEvo.Population(seq, 100)

    new_change = [[10, 1], [20, 2]]
    new = pop.add_sequence(changes=new_change)
    assert (
        new == 100
    ), "new sequence has wrong ID assigned. is {}, should be 100".format(new)
    assert ([10, 1] in pop.changes[new]) & (
        [20, 2] in pop.changes[new]
    ), "changes weren't added correctly. is {}, should be {}".format(
        pop.changes[new], new_change
    )


# consensus_sequence
def test_Population_consensus_sequence_CORRECT():
    seq = SeqSimEvo.Sequence(seq="AAAA")
    pop = SeqSimEvo.Population(seq, 10)
    for _ in range(10):
        pop.add_sequence(changes=[[1, 2]])
    assert (
        str(pop.consensus_sequence()) == "ATAA"
    ), "consensus sequence wrong. is {}, should be {}".format(
        pop.consensus_sequence(), "ATAA"
    )


# copy
def test_Population_copy_CORRECT():
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
def test_Population_delete_sequence_CORRECT():
    seq = SeqSimEvo.Sequence(seq="AAAAAAAAAA")
    pop = SeqSimEvo.Population(seq, 1)

    id1 = pop.add_sequence(changes=[[2, 2]])
    id2 = pop.add_sequence(changes=[[3, 3]])

    pop.delete_sequence(id1)

    assert [3, 3] in pop.get_seq(id1), "{}".format(id2)
    assert pop.n_seq == 2


# get_base
def test_Population_get_base_CORRECT():
    seq = SeqSimEvo.Sequence(seq="AAAAAAAAAA")
    pop = SeqSimEvo.Population(seq, 100)

    assert pop.get_base(0, 1) == 0


# get_sample
def test_Population_get_sample_CORRECT():
    seq = SeqSimEvo.Sequence(seq="AAAAAAAAAA")
    pop = SeqSimEvo.Population(seq, 100)

    sample = pop.get_sample(10)
    assert len(sample) == 10, "sample returns wrong length"
    assert len(set(sample)) == 10, "sample returns double seqIDS"

    sample = pop.get_sample(500)
    assert len(sample) == 100, "sample returns wrong length when oversampling"
    assert len(set(sample)) == 100, "sample returns double seqIDS when oversampling"


# get_seq
def test_Population_get_seq_CORRECT():
    seq = SeqSimEvo.Sequence(seq="AAAAAAAAAA")
    pop = SeqSimEvo.Population(seq, 100)

    id = pop.add_sequence(changes=[[5, 2]])

    assert [5, 2] in pop.get_seq(id)
    assert pop.get_seq(0) is None

    with pytest.raises(Exception) as e_info:
        pop.get_seq(200)


# print_sample
def test_Population_print_sample_correct(capfd):
    seq = SeqSimEvo.Sequence(seq="AAAAAAAAAA")
    pop = SeqSimEvo.Population(seq, 100, changes={1: [[5, 2], [3, 1]]}, changed=[1])

    pop.print_sample([1])
    out, err = capfd.readouterr()

    assert ("A-5-T\t1" in out) and (
        "A-3-G\t1" in out
    ), "printed output:\n{}\nis incomplete. should contain A-5-T\t1 and A-3-G\t1".format(
        out, ""
    )

    with pytest.raises(Exception) as e_info:
        pop.print_sample([200])


# sample_to_string
def test_Population_sample_to_string_correct():
    seq = SeqSimEvo.Sequence(seq="AAAAAAAAAA")
    pop = SeqSimEvo.Population(seq, 100, changes={1: [[5, 2], [3, 1]]}, changed=[1])

    string = pop.sample_to_string([1])
    assert ("A-5-T\t1" in string) and (
        "A-3-G\t1" in string
    ), "returned string:\n{}\nis incomplete. should contain A-5-T\t1 and A-3-G\t1".format(
        string, ""
    )


# stats
def test_Population_stats_correct():
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


# to_fasta
def test_Population_to_fasta_correct():
    seq = SeqSimEvo.Sequence(seq="AAAAAAAAAA")
    pop = SeqSimEvo.Population(seq, 3, changes={1: [[5, 2], [3, 1]]}, changed=[1])

    fasta = ">0\nAAAAAAAAAA\n>1\nAAAGATAAAA\n>2\nAAAAAAAAAA\n"
    assert pop.to_fasta([0, 1, 2]) == fasta, "fasta incorrect"


# Hamming_distance
# def test_Population_hamming_distance_correct():
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


#                           .Simulation.__init__
#                                      .__str__
#                                      .copy
#                                      .get_fitness_effect
#                                      .get_nr_offspring
#                                      .mutate_seq
#                                      .new_generation
#                                      .new_mutations_per_seq

# __init__
"""
def test_Simulation_init_Correct():
    sim = SeqSimEvo.Simulation(
        model="lognormal",
        parameters={"fl": 0.3, "mu": 0.6, "sigma": 0.4},
        mut_rate=0.001,
        subs_matrix=[
            [0, 0.1, 0.6, 1],
            [0.2, 0.2, 0.6, 1],
            [0.4, 0.5, 0.5, 1],
            [0.3, 0.6, 1, 1],
        ],
        seq_len=10,
        basedist=[0, 0.5, 0.8, 1],
        R0=6.5,
        max_pop=10000,
        name="testingSim",
    )

    sim.settings
    assert sim.mut_rate == 0.001
    assert len(sim.sequence) == 10

    assert sim.settings["basedist"] == [0, 0.5, 0.8, 1]
    assert "sequence" not in list(sim.settings.keys())

    assert sim.sequence.sequence.count(0) == 0
    assert sim.gen == 0
    assert sim.average_fitness == 1
    assert sim.effective_pop == 1
    assert sim.n_seq == 1

    assert sim.current_gen.changes == {}
    sim.mutations_per_seq


def test_Simulation_get_fitness_effect():
    seq = SeqSimEvo.Sequence(seq="AAAAAAAAAA")
    sim = SeqSimEvo.Simulation(n_seq_init=3, sequence=seq)

    assert sim.get_fitness_effect(0, 0) == 1

    with pytest.raises(Exception) as e_info:
        sim.get_fitness_effect(50, 0) == 1

    with pytest.raises(Exception) as e_info:
        sim.get_fitness_effect(5, 10) == 1


# copy
def test_Simulation_copy():
    seq = SeqSimEvo.Sequence(seq="AAAAAAAAAA")
    sim = SeqSimEvo.Simulation(n_seq_init=3, sequence=seq)

    sim_copy = sim.copy("sim_test_copy")

    sim.new_generation()

    assert sim.gen != sim_copy.gen


# get_nr_offspring
def test_Simulation_get_nr_offspring():
    seq = SeqSimEvo.Sequence(seq="AAAAAAAAAA")
    sim = SeqSimEvo.Simulation(n_seq_init=3, sequence=seq)
    offspring = sim.get_nr_offspring(0)
    assert offspring >= 0

    offspring, fitness = sim.get_nr_offspring(0, return_fitness=True)
    assert fitness == 1


# mutate_seq
def test_Simulation_mutate_seq():
    seq = SeqSimEvo.Sequence(seq="AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
    sim = SeqSimEvo.Simulation(n_seq_init=3, sequence=seq)
    next_gen = SeqSimEvo.Population(simulation=sim)

    sim.mutations_per_seq = iter([3])
    sim.mutate_seq(next_gen, 0, 0)

    assert next_gen.stats()["total_mutations"] == 3


def test_Simulation_new_generation():
    seq = SeqSimEvo.Sequence.generate_sequence(seq_len=5000)
    sim = SeqSimEvo.Simulation(n_seq_init=10000, sequence=seq, mut_rate=0.001)
    sim.effective_pop = -1
    sim.average_fitness = -1
    sim.n_seq = -1
    old_gen = deepcopy(sim.current_gen.changes)

    sim.new_generation()

    assert sim.gen == 1
    assert (
        sim.current_gen.changes != old_gen
    ), "this test can fail by chance. Be worried if it keeps failing."
    assert sim.effective_pop != -1
    assert sim.average_fitness != -1
    assert sim.n_seq != -1
"""
