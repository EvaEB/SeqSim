"""Sequence core class."""

import random
import numpy as np


class Sequence:
    """Sequence core class.

    Attributes
    ----------
    translation : str
        Order of translation from bases (A, T, G, C) to numbers
    sequence : list
        A string representation of the sequence
    """

    def __init__(self, seq: list):
        """Create a sequence object.

        Parameters
        ----------
        seq : list
            Sequence
        """
        self.translation = "AGTC"

        if isinstance(seq, str):
            seq = self.translate_sequence(seq)

        self.sequence = seq

    @classmethod
    def generate_sequence(cls, seq_len: int, base_dist: list = None):
        """Generate random sequence.

        Parameters
        ----------
        seq_len : int
            Length of sequence.
        base_dist : list
            Base distribution as list in order A, G, T, C. Default: None
        """
        if seq_len < 0:
            raise TypeError("seq_len must be positive integer")
        seq = []
        if base_dist is None:
            base_dist = np.array([0.25, 0.5, 0.75, 1])
        else:
            base_dist = np.array(base_dist)
        for _ in range(seq_len):
            seq.append(4 - sum(base_dist > random.random()))
        return cls(seq)

    def translate_sequence(self, seq: list):
        """Translate a sequences from bases to numbers.

        Parameters
        ----------
        seq : str
            String representation of a sequence.

        Returns
        -------
        list
            Representation of the sequence.
        """
        return [self.translation.index(symbol) for symbol in seq]

    def __getitem__(self, key):
        return self.translation[self.sequence[key]]

    def __str__(self):
        return "".join([self.translation[idx] for idx in self.sequence])

    def __len__(self):
        return len(self.sequence)
