"""Simulation core class."""

import random

from copy import deepcopy
from collections import Counter

import numpy as np
from numpy.typing import ArrayLike

from .sequence import Sequence
from .population import Population
from .simulation_settings import SimulationSettings


class Simulation:
    """Sequence simulation object

    Attributes
    ----------
    sequence : Seq
    settings : SimulationSettings
        Simulation settings.
    fitness_table : ArrayLike
        Table with fitness values for all possible bases in the sequence.
    generation : int
        Current generation.
    average_fitness : int
        Average fitness of the population.
    effective_population : int
        Effective population size of the current population
    n_seq : int
        Total number of sequences in the current population
    """

    def __init__(
        self,
        sequence: Sequence,
        simulation_settings: SimulationSettings,
        fitness_table: ArrayLike = None,
        name: str = "SeqSim",
        **kwargs,
    ):
        """Create a Simulation object.

        Parameters
        ----------
        simulation_settings : SimulationSettings
            Settings for the simulation.
        """
        self.settings = simulation_settings
        self.sequence = sequence
        self.fitness_table = fitness_table
        self.name = name

        if self.fitness_table is None:
            self.fitness_table = self.get_fitness_table()

        self.gen = 0
        self.n_seq = self.settings.n_seq_init

        self.current_population = Population(
            self.sequence, int(self.settings.n_seq_init)
        )

        self._future_mutation_counts = None

    @property
    def effective_population(self):
        """Effective population size."""
        return len(self.current_population)

    @property
    def average_fitness(self):
        """Average fitness of current population."""
        return np.mean(
            [self.get_sequence_fitness(i) for i in range(len(self.current_population))]
        )

    def get_fitness_table(self):
        """Creates a table with random fitness.

        Initiate an array with random fitness values according to the model and
        parameters set in settings.
        """
        if self.settings.model == "neutral":  # neutral model
            fitness = self._get_fitness_neutral()
        elif self.settings.model == "exponential":  # lethals+beneficial+deleterious
            fitness = self._get_fitness_exponential()
        elif self.settings.model == "spikes":
            fitness = self._get_fitness_spikes()
        elif self.settings.model in ["lognormal", "lognormal_ben_only"]:
            fitness = self._get_fitness_lognormal()
        elif self.settings.model == "beta":
            fitness = self._get_fitness_beta()
        elif self.settings.model == "from_data":
            fitness = self._get_fitness_from_data()

        return fitness

    def _get_fitness_neutral(self):
        seq_len = len(self.sequence)
        return np.ones((4, seq_len))

    def _get_fitness_exponential(self):
        seq_len = len(self.sequence)
        fitness = np.zeros((4, seq_len))
        for (idx, base) in enumerate(self.sequence.sequence):
            fitness[base, idx] = 1
        to_fill = np.where(fitness == 0)
        for i, j in zip(*to_fill):
            random_number = random.random()
            if random_number < self.settings.parameters["fl"]:
                fitness[i, j] = 0
            elif (
                random_number
                < self.settings.parameters["fl"] + self.settings.parameters["fn"]
            ):
                fitness[i, j] = 1
            elif (
                random_number
                < self.settings.parameters["fl"]
                + self.settings.parameters["fn"]
                + self.settings.parameters["fb"]
            ):
                fitness[i, j] = 1 + np.random.exponential(
                    self.settings.parameters["lb"]
                )
            else:
                fitness[i, j] = 1 - np.random.exponential(
                    self.settings.parameters["ld"]
                )
            if fitness[i, j] < 0:
                fitness[i, j] = 0
        return fitness

    def _get_fitness_spikes(self):
        seq_len = len(self.sequence)
        fitness = np.zeros((4, seq_len))
        for (idx, base) in enumerate(self.sequence.sequence):
            fitness[base, idx] = 1
        to_fill = np.where(fitness == 0)
        n_spikes = len(self.settings.parameters["loc"])
        for i, j in zip(*to_fill):
            random_number = random.random()
            prob = 0
            for spike in range(n_spikes):
                prob += self.settings.parameters["freq"][spike]
                if random_number < prob:
                    fitness[i, j] = self.settings.parameters["loc"][spike]
                    break
        return fitness

    def _get_fitness_lognormal(self):
        seq_len = len(self.sequence)
        fitness = np.zeros((4, seq_len))
        for (idx, base) in enumerate(self.sequence.sequence):
            fitness[base, idx] = 1
        to_fill = np.where(fitness == 0)
        for i, j in zip(*to_fill):
            random_number = random.random()
            if random_number > self.settings.parameters["fl"]:
                fitness[i, j] = np.random.lognormal(
                    self.settings.parameters["mu"],
                    self.settings.parameters["sigma"],
                )
        if self.settings.model == "lognormal_ben_only":
            fitness[fitness < 1] = 1
        return fitness

    def _get_fitness_beta(self):
        seq_len = len(self.sequence)
        fitness = np.zeros((4, seq_len))
        for (idx, base) in enumerate(self.sequence.sequence):
            fitness[base, idx] = 1
        to_fill = np.where(fitness == 0)
        for i, j in zip(*to_fill):
            random_number = random.random()
            if random_number > self.settings.parameters["fl"]:
                fitness[i, j] = np.random.lognormal(
                    self.settings.parameters["a"], self.settings.parameters["b"]
                )
        return fitness

    def _get_fitness_from_data(self):
        seq_len = len(self.sequence)
        fitness = np.zeros((4, seq_len))
        for (idx, base) in enumerate(self.sequence.sequence):
            fitness[base, idx] = 1
        to_fill = np.where(fitness == 0)
        for i, j in zip(*to_fill):
            index = np.random.randint(len(self.settings.parameters["values"]))
            if self.settings.parameters["SD"][index] == 0:
                new_number = self.settings.parameters["values"][index]
            else:
                new_number = np.random.normal(
                    self.settings.parameters["values"][index],
                    self.settings.parameters["SD"][index],
                )
            if new_number >= 0:
                fitness[i, j] = new_number
            else:
                fitness[i, j] = 0
        return fitness

    def get_fitness_effect(self, location, target_base):
        """Get the fitness effect of a mutation at location to target_base

        Parameters
        ----------
        location : int
        target_base : int

        Returns
        -------
        Fitness effect at this position
        """
        return self.fitness_table[target_base, location]

    def get_sequence_fitness(self, sequence_id):
        """Get fitness for sequence."""
        changes = self.current_population.get_seq(sequence_id)
        fitness = 1

        # calculate the fitness of the sequence
        if changes is not None:
            for pos, base in zip(changes[:, 0], changes[:, 1]):
                fitness *= self.fitness_table[int(base), int(pos)]

        return fitness

    def get_offspring_count(self, fitness):
        """Get number of offsprint of a sequence.

        Get the number of offspring depending on the fitness of that sequence.

        Parameters
        ----------
        fitness : float

        Returns
        -------
        Tuple with the number of offspring and the fitness of the sequence.
        The number of offspring that the sequence produces during one generation.
        """
        basic_reproductive_number = self.settings.basic_reproductive_number

        if self.settings.offspring_distribution == "poisson":
            offspring = np.random.poisson(fitness * basic_reproductive_number)
        elif self.settings.offspring_distribution == "normal":
            offspring = int(
                round(
                    np.random.normal(
                        loc=fitness * basic_reproductive_number,
                        scale=self.settings.offspring_sigma,
                    )
                )
            )
        else:
            raise ValueError(
                f"Offspring distribution {self.settings.offspring_distribution} unknown."
            )
        return offspring

    def mutate_seq(self, pop, seq_id_new, seq_id_old):
        """Mutate a sequence.

        Parameters
        ----------
        pop : Population
            Population the sequence mutates in
        seq_id_new : int
            Sequence ID of the new sequence
        seq_id_old : int
            Sequence ID of the sequence in its ancestor population

        Returns
        -------
        Return True if the sequence did mutate and False if it did not.
        """
        # get the number of mutations that will take place
        mutation_count = self._get_mutation_count()

        if mutation_count <= 0:
            return False

        success_mut = 0
        while success_mut < mutation_count:  # do the mutations one by one
            where = random.randrange(
                0, len(self.sequence)
            )  # draw where the mutation will take place
            base = self.current_population.get_base(seq_id_old, where)
            rand_nr = random.random()  # draw a random nr for base substitution

            to_check = self.settings.subs_matrix[int(base), :]
            # get the cumulative distribution of base substitutions
            new_base = np.where(to_check > rand_nr)[0][0]  # find the new base

            # only accept the mutation if there was an actual change and make
            # sure mutations are accepted in line with the G-A increase
            if base != new_base:
                if (base == 1) and (new_base == 0):  # G-A mutations
                    if (self.settings.ga_increase >= 1) or (
                        random.random() < self.settings["ga_increase"]
                    ):
                        pop.add_change(seq_id_new, where, new_base)
                        success_mut += 1
                elif (base != 1) or (new_base != 0):  # not G-A mutation
                    if (self.settings.ga_increase <= 1) or random.random() < (
                        1.0 / self.settings.ga_increase
                    ):
                        pop.add_change(seq_id_new, where, new_base)
                        success_mut += 1

        return True

    def new_generation(self, new_gen=None):
        """Create a new generation in the simulation.

        Updates state of the simulation.

        Parameters
        ----------
        new_gen : Population
            Population to store the new generation in. If `None` a new
            Population will be created. Default: None.

        Returns
        -------
        The new Population
        """
        self.gen += 1
        if new_gen is None:
            new_gen = Population(self.sequence, n_seq=0)
        all_offspring = []

        fitnesses = [
            self.get_sequence_fitness(i) for i in range(len(self.current_population))
        ]
        weights = [self.get_offspring_count(fitness) for fitness in fitnesses]

        if sum(weights) > self.settings.max_pop:
            # reduce the population randomly to max_pop
            all_offspring = sorted(
                np.random.choice(
                    list(range(self.current_population.n_seq)),
                    size=int(self.settings.max_pop),
                    p=np.array(weights, dtype=float) / sum(weights),
                )
            )
        else:
            all_offspring = [i for i, j in enumerate(weights) for k in range(j)]

        offspring_counter = Counter(all_offspring)
        for seq_id, count in offspring_counter.items():
            changes = self.current_population.get_seq(seq_id)
            for _ in range(count):
                new_seq_id = new_gen.add_sequence(changes=changes)
                self.mutate_seq(new_gen, new_seq_id, seq_id)

        self.current_population = new_gen

        self.n_seq = self.current_population.n_seq

        return self.current_population

    def _get_mutation_count(self):
        """Get a mutation count.

        Uses a pregenerated iterator with mutation counts for efficiency.
        """
        if self._future_mutation_counts is None:
            self._future_mutation_counts = self._get_future_mutation_counts()
        try:
            mutation_count = next(self._future_mutation_counts)
        except StopIteration:
            self._future_mutation_counts = self._get_future_mutation_counts()
            mutation_count = next(self._future_mutation_counts)
        return mutation_count

    def _get_future_mutation_counts(self):
        """Predraw the number of mutations."""
        return iter(
            np.random.binomial(len(self.sequence), self.settings.mutation_rate, 1000)
        )

    def __str__(self):
        string = "sequence simulation\n"
        string += f"MED model\t{self.settings.model}\n"
        for parameter, value in self.settings.parameters.items():
            string += f"{parameter}\t{value}\n"
        string += f"ancestor\t{str(self.sequence)}\n"
        string += f"number of generations\t{str(self.gen)}\n"
        stats = self.current_population.stats()
        for stat, value in stats.items():
            string += f"{stat}\t{value}\n"
        return string

    def copy(self, n_seq=-1, **kwargs):
        """Create a copy of the simulation with settings.

        Parameters
        ----------
        n_seq : int
            Number of sequences to keep. Default: -1.

        Returns
        -------
        Copy of the simulation in its original state.
        """
        if n_seq == -1:  # original state
            return Simulation(
                self.sequence,
                deepcopy(self.settings),
                fitness_table=self.fitness_table,
                **kwargs,
            )
        else:
            simcopy = Simulation(
                self.sequence,
                deepcopy(self.settings),
                fitness_table=self.fitness_table,
                n_seq_init=n_seq,
                **kwargs,
            )
            sample = self.current_population.get_sample(n_seq)
            for i, s in enumerate(sample):
                changes = self.current_population.get_seq(s)
                if changes is not None:
                    for c in changes:
                        simcopy.current_population.add_change(i, c[0], c[1])
            return simcopy
