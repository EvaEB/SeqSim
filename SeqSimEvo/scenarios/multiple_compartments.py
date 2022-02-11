"""multiple_compartments.py: a scenario for SeqSimEvo

Simulation of sequence evolution in multiple compartments, with migration between
the compartments.
"""

import argparse
import ast

import numpy as np

from SeqSimEvo import SimulationSettings, Sequence
from SeqSimEvo.simulation import Simulation


class MultipleCompartments:
    """Simulation of evolution in multiple compartments."""

    def __init__(
        self,
        simulation_settings: SimulationSettings,
        n_compartments: int = 2,
        migration_rate: float = 0.0,
        **kwargs,
    ):
        """Initialize the simulation.

        Parameters
        ----------
        simulation_settings : SimulationSettings
        n_compartments : int
            Number of compartments.
        migration_rate : float | ArrayLike
            Migration rate between the compartments. Default: 0.0.
        names : list
            Names of the compartments. Default: None.
        """
        self.n_compartments = n_compartments

        if isinstance(migration_rate, float):
            self.migration_matrix = (
                np.ones((n_compartments, n_compartments)) * migration_rate
            )
        else:
            self.migration_matrix = migration_rate

        # set up sims for all compartments
        sequence = Sequence.generate_sequence(
            simulation_settings.seq_len, simulation_settings.base_dist
        )
        self.sims = [
            Simulation(sequence, simulation_settings, name=str(i))
            for i in range(self.n_compartments)
        ]

    def new_generation(self):
        """preform migration and generate a new generation for all compartments"""
        # do migration
        # get emigrants from each compartment
        migrants = [[] for i in range(self.n_compartments)]
        for comp1 in range(self.n_compartments):
            for comp2 in range(self.n_compartments):
                if comp1 != comp2:
                    # select migrants
                    mig_rate = self.migration_matrix[comp1, comp2]
                    n_to_migrate = mig_rate * self.sims[comp1].current_population.n_seq
                    n_to_migrate = int(n_to_migrate) + (
                        1
                        if np.random.random() < (n_to_migrate - int(n_to_migrate))
                        else 0
                    )
                    to_migrate = self.sims[comp1].current_population.get_sample(
                        n_to_migrate
                    )

                    # remember changes in these migrants for later adding, remove from the current compartment
                    for seq in to_migrate:
                        migrants[comp2].append(
                            self.sims[comp1].current_population.get_seq(seq)
                        )
                        self.sims[comp1].current_population.delete_sequence(seq)

        # add immigrants to each compartment
        for comp in range(self.n_compartments):
            for change in migrants[comp]:
                self.sims[comp].current_population.add_sequence(changes=change)

            # let this compartment generate a new generation
            self.sims[comp].new_generation()
        if sum([i.n_seq for i in self.sims]) == 0:
            dieout = True
            for i in self.sims:
                i.n_seq = i.settings["n_seq_init"]
            self.new_generation()

    def __str__(self):
        out = ""
        for sim in self.sims:
            out += f"compartment {sim.name}\t{sim.current_population.n_seq} seqs\n"
        return out


def run(scenario_settings, preset):
    kw_settings = {}
    for i in scenario_settings:
        if i in [
            "n_compartments",
            "diverse_index",
            "names",
            "n_seq_init",
            "migration_rate",
            "mut_rate",
            "R0",
            "max_pop",
        ]:
            kw_settings[i] = scenario_settings[i]
    simulation_settings = SimulationSettings.from_preset(preset)
    sim = MultipleCompartments(simulation_settings, **kw_settings)

    if type(scenario_settings["sampling_amount"]) is int:
        n_seq = [scenario_settings["sampling_amount"]] * scenario_settings[
            "n_compartments"
        ]
    else:
        n_seq = scenario_settings["sampling_amount"]

    fasta = ""

    for i in range(scenario_settings["n_gen"]):
        print(i)
        sim.new_generation()
        if i + 1 in scenario_settings["sampling_times"]:
            for j, s in enumerate(sim.sims):
                fasta += s.current_population.to_fasta(
                    n_seq=n_seq[j],
                    description=f"-{s.settings.name}-gen{i+1}",
                )

    return fasta


def main():
    parser = argparse.ArgumentParser(
        description="simulation of evolution in multiple compartments"
    )

    parser.add_argument(
        "-o", type=str, default="HIV", help="organism settings to use, defaults to HIV"
    )
    parser.add_argument(
        "-nc",
        type=int,
        default=2,
        help="number of compartments to simulate. Defaults to 2",
    )
    parser.add_argument(
        "-ng",
        type=int,
        default=100,
        help="number of generations to simulate. Defaults to 100",
    )
    parser.add_argument(
        "-d",
        type=float,
        default=0,
        help="how different the fitness tables are for each of \
                                  the compartments (0 - all the same, parallel \
                                  evolution, 1 - all different, complete divergent \
                                  evolution). Defaults to 0",
    )
    parser.add_argument(
        "-mig",
        default="[[0,0],[0,0]]",
        type=str,
        help="migration rate matrix, \
                                  or single value (use migration rate between all populations),\
                                  defaults to no migration",
    )
    parser.add_argument(
        "-names", nargs="+", default=None, help="names of the compartments"
    )
    parser.add_argument(
        "-ninit",
        nargs="+",
        default=None,
        type=int,
        help="number of initial sequences per compartment, defaults to one per compartment",
    )
    parser.add_argument(
        "-mut",
        nargs="+",
        default=None,
        type=float,
        help="mutation rates per compartment, defaults to the value in the organism settings",
    )
    parser.add_argument(
        "-R0", nargs="+", default=None, type=float, help="R0 per compartment"
    )
    parser.add_argument(
        "-maxpop",
        nargs="+",
        default=None,
        type=int,
        help="maximum population per compartment",
    )
    parser.add_argument(
        "-st",
        nargs="+",
        default=None,
        type=int,
        help="sampling times per compartment,defaults to last",
    )
    parser.add_argument(
        "-sa",
        nargs="+",
        default=10,
        type=int,
        help="sampling amount per compartment,defaults to 10",
    )

    args = parser.parse_args()

    settings = {}
    settings["n_compartments"] = args.nc
    settings["diverse_index"] = args.d
    settings["n_gen"] = args.ng

    if "[" not in args.mig:
        settings["migration"] = np.ones([len(args.ninit)] * 2) * float(args.mig[0])
    else:
        settings["migration"] = ast.literal_eval(args.mig)

    settings["sampling_amount"] = args.sa

    if args.st is None:
        settings["sampling_times"] = [args.ng]
    else:
        settings["sampling_times"] = args.st

    if args.names is not None:
        settings["names"] = args.names

    if args.ninit is not None:
        settings["n_seq_init"] = args.ninit

    if args.mut is not None:
        settings["mut_rate"] = args.mut

    if args.R0 is not None:
        settings["R0"] = args.R0

    if args.maxpop is not None:
        settings["max_pop"] = args.maxpop

    print(run(settings, args.o))


if __name__ == "__main__":
    main()
