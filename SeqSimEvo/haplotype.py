"""Haplotype."""

from copy import deepcopy, copy

import numpy as np
from numpy.typing import ArrayLike

from .sequence import Sequence


class Haplotype:
    """Haplotype core class."""

    def __init__(self, sequence: Sequence, changes: ArrayLike = None):
        """Initialize haplotype."""
        self.sequence = sequence
        self.changes = changes

    @property
    def mutation_count(self):
        """Get number of mutations."""
        return len(self.changes)

    def get_base(self, pos: int):
        """Get specific base for a sequence.

        Parameters
        ----------
        pos : int
        """
        if self.changes is not None and pos in self.changes[:, 0]:
            return self.changes[self.changes[:, 0] == pos, 1]
        return self.sequence.sequence[pos]

    def set_base(self, pos: int, base: int):
        """Set specific base for a sequence.

        Parameters
        ----------
        pos : int
        base : int
        """
        if self.changes is None:
            self.changes = np.array([[pos, base]])
            return

        if pos in self.changes[:, 0]:
            self.changes[self.changes[:, 0] == pos, 1] = base
            return

        self.changes = np.vstack((self.changes, [pos, base]))

    def get_sequence(self):
        """Get sequence with mutations."""
        if self.changes is None:
            return self.sequence

        sequence = deepcopy(self.sequence)
        for pos, base in self.changes:
            sequence.sequence[pos] = base
        return sequence

    def copy(self):
        """Copy haplotype."""
        return Haplotype(self.sequence, copy(self.changes))

    def __getitem__(self, key):
        return self.get_base(key)
