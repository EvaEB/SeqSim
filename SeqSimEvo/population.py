"""Population core class."""

import random

from copy import deepcopy
from collections import Counter
from typing import List

import numpy as np
import scipy as sp

from numpy.typing import NDArray

from .sequence import Sequence
from .haplotype import Haplotype


class Population:
    """Population core class.

    Attributes
    ----------
    sequence : Sequence
        Reference sequence
    n_seq : int
        Number of sequences in the population
    haplotypes : ndarray
        Maps sequence ids to an array that stores all mutations
    changed : set
        Sequence ids that contain mutations
    """

    def __init__(self, sequence: Sequence, n_seq: int, haplotypes: NDArray = None):
        """Initialize the Population

        Parameters
        ----------
        sequence : Sequence
        n_seq : int
            Population size.
        """
        self.sequence = sequence
        self.n_seq = n_seq

        if haplotypes is None:
            haplotypes = np.full(n_seq, Haplotype(self.sequence))

        self.haplotypes = haplotypes

    @property
    def changed(self):
        """Get the mutated sequences."""
        return [
            seq_id
            for seq_id, haplotype in enumerate(self.haplotypes)
            if haplotype.changes is not None
        ]

    @classmethod
    def merge(cls, *populations):
        """Merge a list of populations."""
        sequence = populations[0].sequence
        if not all(sequence == population.sequence for population in populations):
            raise ValueError("Can only merge Populations with same reference sequence.")

        n_seq = sum([len(pop) for pop in populations])
        new_population = cls(sequence, n_seq)
        seq_id = 0
        for population in populations:
            new_population.haplotypes[
                seq_id : seq_id + len(population)
            ] = population.haplotypes[:]
            seq_id += len(population)
        return new_population

    def __add__(self, other):
        if self.sequence != other.sequence:
            raise ValueError("Can only add Populations with same reference sequence.")
        return Population.merge(self, other)

    def split(self, ratios: list):
        """Split Population into segments."""
        ratios = np.array(ratios) / sum(ratios)
        sizes = [int(ratio * self.n_seq) for ratio in ratios]
        if sum(sizes) > self.n_seq:
            sizes[-1] -= 1
        if sum(sizes) < self.n_seq:
            sizes[-1] += 1
        seq_id = 0
        populations = []
        for size in sizes:
            population = Population(
                self.sequence, size, haplotypes=self.haplotypes[seq_id : seq_id + size]
            )
            seq_id += size
            populations.append(population)
        return populations

    def get_haplotype(self, seq_id: int) -> Haplotype:
        """Get changes for sequence."""
        if seq_id < 0 or seq_id >= self.n_seq:
            raise IndexError(f"Sequence id {seq_id} out of range.")
        return self.haplotypes[seq_id]

    def get_base(self, seq_id: int, pos: int):
        """Get specific base for a sequence.

        Parameters
        ----------
        seq_id : int
        pos : int
        """
        return self.haplotypes[seq_id].get_base(pos)

    def get_sample(self, sample_size: int) -> List[int]:
        """Get a random sample of the population.

        If the sample size is larger than the population, the whole population
        is returned.

        Paramters
        ---------
        sample_size : int

        Returns
        -------
        list
            List of randomly selected sample sequence ids.
        """
        if sample_size >= self.n_seq:
            return list(range(self.n_seq))
        return np.random.choice(self.n_seq, size=int(sample_size), replace=False)

    def set_haplotype(self, seq_id: int, haplotype: Haplotype):
        """Set haplotype for seq_id."""
        self.haplotypes[seq_id] = haplotype

    def add_change(self, seq_id: int, pos: int, target: int):
        """Add a change to an existing sequence.

        Parameters
        ----------
        seq_id : int
        pos : int
            Position of the change.
        target : int
            Target base.
        """
        if pos > len(self.sequence):
            raise IndexError(f"{pos=} is outsite of sequence length.")

        if seq_id > self.n_seq:
            raise IndexError(f"{seq_id=} is not in the population.")

        haplotype = self.haplotypes[seq_id].copy()
        haplotype.set_base(pos, target)
        self.haplotypes[seq_id] = haplotype

    def stats(self) -> dict:
        """Get stats for the population.

        Returns
        -------
        dict
            * `n_seq`: the total number of sequences in the current generation
            * `unmutated`: the number of unmutated sequences
            * `total_mutations`: the number of mutations in total
            * `unique_mutations`: the length of the set of all mutations
            * `majority_mutations`: the number of mutations that reached majority
            * `max_fraction`: the highest fraction reached by a mutation
            * `GA_rate`: the fraction of mutations that are G-to-A
        """
        stats = {}
        stats["n_seq"] = self.n_seq
        stats["unmutated"] = self.n_seq - len(self.changed)
        if len(self.changed) > 0:
            all_mutations = np.vstack(
                list(haplotype.changes for haplotype in self.haplotypes[self.changed])
            )
        else:
            all_mutations = []
        stats["total_mutations"] = len(all_mutations)
        all_mutations = [tuple(row) for row in all_mutations]
        stats["unique_mutations"] = len(set(all_mutations))

        if len(all_mutations) > 0:
            mut_counts = np.array(list(Counter(all_mutations).values()))
        else:
            mut_counts = []
        if len(mut_counts) > 0:
            stats["majority_mutations"] = sum(mut_counts > (stats["n_seq"] / 2.0))
            stats["max_fraction"] = max(mut_counts / float(stats["n_seq"]))
        else:
            stats["majority_mutations"] = 0
            stats["max_fraction"] = 0
        ga = 0
        for i in all_mutations:
            if self.sequence[i[0]] == 1 and i[1] == 0:
                ga += 1.0
        try:
            stats["GA_rate"] = ga / len(all_mutations)
        except ZeroDivisionError:
            stats["GA_rate"] = None
        return stats

    def to_fasta(
        self,
        seq_ids: List[int] = None,
        n_seq: int = None,
        description: str = "",
    ) -> str:
        """Convert (part of) the population to fasta-format.

        Without any arguments, all sequences in the population will be returned.

        Parameters
        ----------
        seq_ids : list
            List of sequence ids to include.
        n_seq : int
            Number of sequences to convert to fasta-format. Drawn at random.
            Will be ignored if `seq_ids` is not None.
        description : str
            Description of the sequences will be added after the sequence id
            in the name of each sequence.

        Returns
        -------
        str
            Selected sequences in fasta-format.
        """
        if seq_ids is None:
            seq_ids = []
        string = ""
        if len(seq_ids) == 0:
            if n_seq is None or n_seq > self.n_seq:
                n_seq = self.n_seq
            seq_ids = random.sample(list(range(self.n_seq)), n_seq)

        for seq_id in seq_ids:
            string += ">" + str(seq_id) + "" + str(description) + "\n"
            string += str(self.get_haplotype(seq_id).get_sequence()) + "\n"
        return string

    def consensus_sequence(self) -> Sequence:
        """Compute consensus sequence of the population."""
        seq = deepcopy(self.sequence)
        changes = np.vstack(
            list(
                haplotype.changes
                for haplotype in self.haplotypes
                if haplotype.changes is not None
            )
        )
        mutations = [tuple(mut) for mut in changes]

        mutation_abundances = Counter(mutations)
        for mut, count in mutation_abundances.items():
            if count >= self.n_seq / 2.0:
                seq.sequence[int(mut[0])] = int(mut[1])
        return seq

    def hamming_distance(self, sample, simulation_settings, action="mean"):
        """Calculate the inter-sequence hamming distances in a sample.

        if action is 'mean', return the mean hamming distance,
        if action is 'Poisson_fit', return the poisson fit for the time since
            infection from the distribution of hamming distances as presented
            in Lee et al, 2010
        """
        HDs = []
        for i in sample:
            if i in self.changed:
                changed1 = [str(k) for k in self.haplotypes[i].changes]
            else:
                changed1 = []
            for j in sample:
                if i != j:
                    if j in self.changed:
                        changed2 = [str(k) for k in self.haplotypes[j].changes]
                    else:
                        changed2 = []
                    HDs.append(len(set(list(changed1)) ^ set(list(changed2))))
        if action == "mean":
            return np.mean(HDs)
        elif action == "Poisson_fit":
            poiss = np.mean(HDs) / (
                2 * simulation_settings["mut_rate"] * simulation_settings["seq_len"]
            )
            exp = sp.stats.poisson.pmf(list(range(max(HDs) + 1)), np.mean(HDs)) * len(
                HDs
            )
            obs = np.histogram(HDs, bins=list(range(0, max(HDs) + 2)))[0]
            pval = sp.stats.chisquare(obs, exp, ddof=len(exp) - 1 - len(sample)).pvalue
            if np.isnan(pval) or pval > 0.05:
                return poiss
            else:
                return np.nan

    def print_sample(self, seq_ids: list):
        """Print a summary of the mutations that have occurred

        Prints a summary of the mutations that have occurred in all seq_ids in
        the format: #mutID (from-pos-to)\tsequence\tpatient\n

        Parameters
        ----------
        seq_ids : list
            List of sequences to print a summary of.
        """
        print(self.get_summary(seq_ids))

    def get_summary(self, seq_ids: list):
        """Create a summary of the mutations that have occurred.

        Parameters
        ----------
        seq_ids : list
            List of sequences to print a summary of.

        Returns
        -------
        string
            Summary of the mutations that have occurred in `seq_ids`, in format
            `#mutID (from-pos-to)\\tsequence\\tpatient`
        """
        string = "#mutID (from-pos-to)\tsequence\tpatient\n"
        changed = self.changed
        for i in sorted(seq_ids):
            if i in changed:
                for pos, base_id in self.haplotypes[i].changes:
                    base = self.sequence.translation[base_id]
                    string += f"{self.sequence[pos]}-{pos}-{base}\t{i}\n"
        return string

    def copy(self):
        """Create a deep copy of the population"""
        return Population(
            sequence=self.sequence,
            n_seq=self.n_seq,
            haplotypes=deepcopy(self.haplotypes),
        )

    def __len__(self):
        return self.n_seq

    def __str__(self):
        return self.get_summary(range(self.n_seq))
