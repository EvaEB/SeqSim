"""Simulation settings module."""

import os
from dataclasses import dataclass

import yaml
from pkg_resources import resource_filename

import numpy as np
from numpy.typing import ArrayLike


@dataclass
class SimulationSettings:
    """Simulation settings."""

    model: str
    offspring_distribution: str

    mutation_rate: float
    subs_matrix: ArrayLike

    seq_len: int
    base_dist: list

    n_seq_init: int
    basic_reproductive_number: float
    ga_increase: int
    max_pop: int
    name: str

    parameters: dict = None
    offspring_sigma: float = 0
    dilution: float = -1.0

    def __post_init__(self):
        if isinstance(self.subs_matrix, list):
            self.subs_matrix = np.array(self.subs_matrix)

    @classmethod
    def from_preset(cls, name: str):
        """Create settings from preset."""
        preset_path = resource_filename("SeqSimEvo", "simulation_settings/")
        try:
            with open(os.path.join(preset_path, name), encoding="utf8") as preset_file:
                settings = yaml.safe_load(preset_file)
        except FileNotFoundError as err:
            raise FileNotFoundError(f"Preset {name} does not exist.") from err
        return cls.from_dict(settings)

    @classmethod
    def from_dict(cls, settings: dict):
        """Create settings from dict."""
        return cls(**settings)
