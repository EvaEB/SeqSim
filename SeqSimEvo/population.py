"""Population core class."""

import random

from copy import deepcopy, copy
from collections import Counter

import numpy as np
import scipy as sp


class Population:
    """Population core class.

    Attributes
    ----------
    sequence : Sequence
        Reference sequence
    changed : set
        Sequence ids that contain mutations
    changes : dict
        Maps sequence ids to an array that stores all mutations
    n_seq : int
        Number of sequences in the population
    """

    def __init__(self, sequence, n_seq, changes=None, changed=None, **kwargs):
        """Initialize the Population

        Arguments:
            changes (dict): for pre-existing changes. Numpy array per sequence with
                the position and new base per mutation
            changed (set): for pre-existing changes: the sequence IDs that countain a mutation
            n_seq (int): number of sequences in the population (default 1)
        """
        self.sequence = sequence
        if changes is None:
            self.changed = set([])
            self.changes = {}
        else:
            self.changes = changes
            for i in self.changes:
                self.changes[i] = np.array(self.changes[i])
            self.changed = set(changed)
        self.n_seq = n_seq

    def copy(self):
        """Create a deep copy of the population"""
        return Population(
            sequence=self.sequence,
            changes=deepcopy(self.changes),
            changed=deepcopy(self.changed),
            n_seq=self.n_seq,
        )

    def print_sample(self, seq_ids):
        """Print a summary of the mutations that have occurred

        Prints a summary of the mutations that have occurred in all seq_ids in
        the format: #mutID (from-pos-to)\tsequence\tpatient\n

        **Arguments**:
        * `seq_ids` (list): the ids of the sequences to print a summary of

        **Returns**:
        Nothing. Output is printed to stdout.
        """
        if any(np.array(seq_ids) > self.n_seq):
            raise IndexError("seqID out of range")
        string = "#mutID (from-pos-to)\tsequence\tpatient\n"
        for i in range(self.n_seq):
            if i in self.changed and i in seq_ids:
                for j in self.changes[i]:
                    pos = j[0]
                    string += "{orig}-{pos}-{to}\t{seq}\t\n".format(
                        orig=self.sequence[pos],
                        pos=pos,
                        to=self.sequence.translation[j[1]],
                        seq=i,
                    )
        print(string)

    def sample_to_string(self, seq_ids):
        """create a summary of the mutations that have occured

        **Arguments**:
        * `seq_ids` (list): the ids of the sequences to print a summary of

        **Returns**:
        * `string`: summary of the mutations that have occured in the seq_ids,
        in the format "#mutID (from-pos-to)\tsequence\tpatient"

        """
        string = "#mutID (from-pos-to)\tsequence\tpatient\n"
        for i in range(self.n_seq):
            if i in self.changed and i in seq_ids:
                for j in self.changes[i]:
                    pos = j[0]
                    string += "{orig}-{pos}-{to}\t{seq}\t\n".format(
                        orig=self.sequence[pos],
                        pos=pos,
                        to=self.sequence.translation[j[1]],
                        seq=i,
                    )
        return string

    def get_sample(self, sample_size):
        """get a random sample from the population
        If the sample size is larger than the population, the whole population is
        returned

        **Arguments**:
        * `sample_size` (int): the size of the sample

        **Returns**:
            a list of sequence IDs randomly sampled from the population
        """
        try:
            return np.random.choice(self.n_seq, size=int(sample_size), replace=False)
        except ValueError:
            return list(range(self.n_seq))

    def delete_sequence(self, ID):
        """delete sequence from the population.
        Sequence IDs will be reassigned to fit within the new number of sequences.

        **Arguments**:
        * `ID` (int): the sequence ID of the sequence to remove

        **Returns**:
        nothing. Sequence is removed in-place.
        """
        self.n_seq -= 1
        if self.get_seq(ID) is not None:
            self.changed.remove(ID)
            del self.changes[ID]

        if self.get_seq(self.n_seq) is not None:
            self.changed.remove(self.n_seq)
            self.changed.add(ID)

            self.changes[ID] = self.changes[self.n_seq]
            del self.changes[self.n_seq]

    def add_sequence(self, changes=None):
        """add a sequence to the population

        add a sequence, optionally with certain changes (as a list of position, new), to the population

        **Arguments**:
        * `changes`: the changes present in this sequence

        **Returns**:
        the sequence ID of the newly added sequence
        """
        self.n_seq += 1
        if changes is not None:
            for i in changes:
                self.add_change(self.n_seq - 1, i[0], i[1])

        return self.n_seq - 1

    def add_change(self, seq_id, pos, target):
        """add a change to an existing sequence

        **Arguments**:
        * `seq_id` (int): the sequence ID to add the change to
        * `pos` (int): the position of the change
        * `target` (int): the new base at the changed position

        **Returns**:
            Nothing. Population is changed in-place.
        """
        if pos > len(self.sequence):
            raise IndexError("Pos {} outside sequence length".format(pos))

        if seq_id > self.n_seq:
            raise IndexError(
                "SeqID {} outside pop size {} {}".format(seq_id, self.n_seq, self)
            )

        if seq_id in self.changed:
            # add to existing changes list
            if pos in self.changes[seq_id][:, 0]:
                self.changes[seq_id][self.changes[seq_id][:, 0] == pos, 1] = target
            else:
                self.changes[seq_id] = np.vstack((self.changes[seq_id], [pos, target]))
        else:
            # add a new changed sequence
            self.changed.add(seq_id)
            self.changes[seq_id] = np.array([[pos, target]])

    def get_base(self, seq_id, pos):
        """get the current base at position pos in sequence with id seq_id"""
        if seq_id in self.changed:
            if pos in self.changes[seq_id][:, 0]:
                return self.changes[seq_id][self.changes[seq_id][:, 0] == pos, 1]
        return self.sequence.sequence[pos]

    def stats(self):
        """return a dict of stats about the population

        **Keys in the returned dict:**
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
        stats["unmutated"] = self.n_seq - len(self.changes)
        if len(self.changed) > 0:
            all_mutations = np.vstack(list(self.changes.values()))
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
        GA = 0
        for i in all_mutations:
            if self.sequence[i[0]] == 1 and i[1] == 0:
                GA += 1.0
        try:
            stats["GA_rate"] = GA / len(all_mutations)
        except ZeroDivisionError:
            stats["GA_rate"] = None
        return stats

    def to_fasta(self, seq_ids=[], n_seq=None, description="", progress=False):
        """
        convert (part of) the population to fasta-format.

        Without any arguments, all sequences in the population will be returned.

        **Optional arguments**:
        * `seq_ids` (list): list of sequence IDs to convert to fasta-format
        * `n_seq`: number of sequences to convert to fasta-format (random draw
        from population). This number will be ignored if seq_ids is given
        * `description` (str): description of the sequences, will be added after
        the sequenceID in the header of each sequence
        * `progress` (Bool): display a progress bar?

        **Returns**:
            str: the selected sequences in fasta-format
        """
        string = ""
        if len(seq_ids) == 0:
            if n_seq is None or n_seq > self.n_seq:
                n_seq = self.n_seq

            seq_ids = random.sample(list(range(self.n_seq)), n_seq)
        for i in range(len(seq_ids)):
            seqID = seq_ids[i]
            string += ">" + str(seqID) + "" + str(description) + "\n"
            changed_here = self.get_seq(seqID)
            seq = deepcopy(self.sequence)
            if changed_here is not None:
                for i in changed_here:
                    seq.sequence[int(i[0])] = int(i[1])
            string += str(seq) + "\n"
        return string

    def consensus_sequence(self):
        """return the consensus sequence of the population"""
        seq = deepcopy(self.sequence)
        all_mutations = np.vstack(list(self.changes.values()))
        all_mutations = [tuple(row) for row in all_mutations]

        mutations = Counter(all_mutations)
        for mut in list(mutations.keys()):
            if mutations[mut] >= self.n_seq / 2.0:
                seq.sequence[int(mut[0])] = int(mut[1])
        return seq

    def get_seq(self, sequence_id):
        """get the changes in the sequence with id sequence_id

        **Raises**:
            `IndexError`: when sequence_id is out of bounds
        """
        if sequence_id > self.n_seq:
            raise IndexError("sequence_id is out of bounds")
        elif sequence_id in self.changed:
            return self.changes[sequence_id]
        else:
            return None

    def hamming_distance(self, sample, simulation_settings, action="mean"):
        """
        calculate the inter-sequence hamming distances in a sample.
        if action is 'mean', return the mean hamming distance,
        if action is 'Poisson_fit', return the poisson fit for the time since
            infection from the distribution of hamming distances as presented
            in Lee et al, 2010
        """
        HDs = []
        for i in sample:
            if i in self.changed:
                changed1 = [str(k) for k in self.changes[i]]
            else:
                changed1 = []
            for j in sample:
                if i != j:
                    if j in self.changed:
                        changed2 = [str(k) for k in self.changes[j]]
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

    def __len__(self):
        return self.n_seq

    def __str__(self):
        string = "#mutID (from-pos-to)\tsequence\tpatient\n"
        for i in range(self.n_seq):
            if i in self.changed:
                for j in self.changes[i]:
                    pos = j[0]
                    string += "{orig}-{pos}-{to}\t{seq}\t\n".format(
                        orig=self.sequence[pos],
                        pos=pos,
                        to=self.sequence.translation[j[1]],
                        seq=i,
                    )
        return string
