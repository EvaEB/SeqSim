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
    fitnesses : list
        Fitnesses of the current population.
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
        **_kwargs,
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

        self.current_population = Population(
            self.sequence, int(self.settings.n_seq_init)
        )

        self._future_mutation_counts = None

    @property
    def effective_population(self):
        """Effective population size."""
        return len(self.current_population)

    @property
    def fitnesses(self):
        """Fitnesses of the current population."""
        return [
            self.get_sequence_fitness(i) for i in range(len(self.current_population))
        ]

    @property
    def average_fitness(self):
        """Average fitness of current population."""
        return np.mean(self.fitnesses)

    @property
    def n_seq(self):
        """Current number of sequences in population."""
        return len(self.current_population)

    @classmethod
    def create_fitness_table(cls, sequence: Sequence, settings: SimulationSettings):
        """Create a fitness table."""
        return cls(sequence, settings).fitness_table

    def get_fitness_table(self):
        """Create a table with random fitness.

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

    def get_sequence_fitness(self, seq_id):
        """Get fitness for sequence."""
        haplotype = self.current_population.get_haplotype(seq_id)
        fitness = 1

        # calculate the fitness of the sequence
        changes = haplotype.changes
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

    def mutate_sequence(self, population: Population, seq_id: int):
        """Mutate a sequence.

        Parameters
        ----------
        pop : Population
            Population the sequence mutates in
        seq_id : int
            Sequence id to mutate.
        """
        mutation_count = self._get_mutation_count()

        if mutation_count <= 0:
            return False

        haplotype = population.get_haplotype(seq_id)
        new_haplotype = haplotype.copy()

        success_mut = 0
        while success_mut < mutation_count:  # do the mutations one by one
            pos = random.randrange(0, len(self.sequence))
            base = haplotype.get_base(pos)

            substitution_vector = self.settings.subs_matrix[int(base), :]
            new_base = np.where(substitution_vector > random.random())[0][0]

            if base == new_base:
                continue

            if (base == 1) and (new_base == 0):  # G-A mutations
                if (self.settings.ga_increase >= 1) or (
                    random.random() < self.settings["ga_increase"]
                ):
                    new_haplotype.set_base(pos, new_base)
                    success_mut += 1
            elif (base != 1) or (new_base != 0):  # not G-A mutation
                if (self.settings.ga_increase <= 1) or random.random() < (
                    1.0 / self.settings.ga_increase
                ):
                    new_haplotype.set_base(pos, new_base)
                    success_mut += 1

        population.set_haplotype(seq_id, new_haplotype)

        return True

    def generate_offspring(self):
        """Generate offspring counts for all sequences."""
        return [self.get_offspring_count(fitness) for fitness in self.fitnesses]

    def subsample_population(self, weights: list, factor: int = 1, sort: bool = True):
        """Subsample the amplified population."""
        pop_size = sum(weights)
        target_size = self.settings.max_pop

        if self.settings.dilution > 0:
            target_size = min(pop_size * self.settings.dilution, self.settings.max_pop)

        if pop_size <= target_size:
            return [i for i, j in enumerate(weights) for k in range(j)]

        sample = np.random.choice(
            list(range(self.current_population.n_seq)),
            size=int(target_size * factor),
            p=np.array(weights, dtype=float) / sum(weights),
        )

        if sort:
            sample = sorted(sample)

        return sample

    def population_from_offspring(self, offspring):
        """Generate population from offspring."""
        new_population = Population(self.sequence, n_seq=len(offspring))
        offspring_counter = Counter(offspring)
        new_seq_id = 0
        for seq_id, count in offspring_counter.items():
            haplotype = self.current_population.get_haplotype(seq_id)
            for _ in range(count):
                new_population.set_haplotype(new_seq_id, haplotype)
                self.mutate_sequence(new_population, new_seq_id)
                new_seq_id += 1
        return new_population

    def new_generation(self):
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
        weights = self.generate_offspring()
        offspring = self.subsample_population(weights)

        self.gen += 1
        self.current_population = self.population_from_offspring(offspring)

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
        new_simulation = Simulation(
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
                    new_simulation.current_population.add_change(i, c[0], c[1])
        return new_simulation
