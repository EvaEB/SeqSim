"""Simulation module tests."""

from SeqSimEvo.simulation import SimulationSettings
import pytest
from copy import deepcopy

import SeqSimEvo


#                           .Simulation.__init__
#                                      .__str__
#                                      .copy
#                                      .get_fitness_effect
#                                      .get_nr_offspring
#                                      .mutate_seq
#                                      .new_generation
#                                      .new_mutations_per_seq

# __init__
def test_simulation_init():
    settings = SeqSimEvo.SimulationSettings.from_preset("HIV")
    seq = SeqSimEvo.Sequence.generate_sequence(100)
    sim = SeqSimEvo.Simulation(seq, settings)

    assert hasattr(sim, "settings")


def test_simulation_get_fitness_effect():
    settings = SeqSimEvo.SimulationSettings.from_preset("HIV")
    settings.n_seq_init = 3
    seq = SeqSimEvo.Sequence(seq="AAAAAAAAAA")
    sim = SeqSimEvo.Simulation(seq, settings)

    assert sim.get_fitness_effect(0, 0) == 1

    with pytest.raises(Exception) as e_info:
        sim.get_fitness_effect(50, 0) == 1

    with pytest.raises(Exception) as e_info:
        sim.get_fitness_effect(5, 10) == 1


# copy
def test_simulation_copy():
    settings = SeqSimEvo.SimulationSettings.from_preset("HIV")
    settings.n_seq_init = 3
    seq = SeqSimEvo.Sequence(seq="AAAAAAAAAA")
    sim = SeqSimEvo.Simulation(seq, settings)

    sim_copy = sim.copy("sim_test_copy")

    sim.new_generation()

    assert sim.gen != sim_copy.gen


# get_nr_offspring
def test_simulation_get_nr_offspring():
    settings = SeqSimEvo.SimulationSettings.from_preset("HIV")
    settings.n_seq_init = 3
    seq = SeqSimEvo.Sequence(seq="AAAAAAAAAA")
    sim = SeqSimEvo.Simulation(seq, settings)

    fitness = sim.get_sequence_fitness(0)
    assert fitness == 1
    offspring = sim.get_offspring_count(fitness)
    assert offspring >= 0


# mutate_seq
def test_simulation_mutate_seq():
    settings = SeqSimEvo.SimulationSettings.from_preset("HIV")
    settings.n_seq_init = 3
    seq = SeqSimEvo.Sequence(seq="AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
    sim = SeqSimEvo.Simulation(seq, settings)
    next_gen = SeqSimEvo.Population(seq, settings.n_seq_init)

    sim._future_mutation_counts = iter([3])
    sim.mutate_seq(next_gen, 0, 0)

    assert next_gen.stats()["total_mutations"] == 3


def test_simulation_new_generation():
    settings = SeqSimEvo.SimulationSettings.from_preset("HIV")
    settings.n_seq_init = 10000
    settings.mutation_rate = 1e-3
    seq = SeqSimEvo.Sequence.generate_sequence(seq_len=5000)
    sim = SeqSimEvo.Simulation(seq, settings)
    sim.n_seq = -1
    old_gen = deepcopy(sim.current_population.changes)

    sim.new_generation()

    assert sim.gen == 1
    assert (
        sim.current_population.changes != old_gen
    ), "this test can fail by chance. Be worried if it keeps failing."
    assert sim.n_seq != -1
