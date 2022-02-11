"""Virus Passaging experiment.

Recreates a virus passaging experiment, where virus is allowed target multiply on a
host for a certain number of generations, after which a proportion is
transferred target new hosts.
"""

import argparse
import ast
from copy import deepcopy as copy
from itertools import permutations

import numpy as np
from tqdm import tqdm

from SeqSimEvo import Sequence, Simulation, SimulationSettings


class Passaging:
    """Passaging experiment."""

    def __init__(self, simulation_settings, passaging_settings):
        self.settings = copy(passaging_settings)
        sequence = Sequence.generate_sequence(simulation_settings.n_seq_init)
        self.sims = []
        self.n_pop = len(self.settings["pop_size"])

        for idx, n_pop in enumerate(self.settings["pop_size"]):
            n_pop_settings = copy(simulation_settings)
            n_pop_settings.name = str(idx)
            n_pop_settings.n_seq_init = n_pop
            self.sims.append(Simulation(copy(sequence), n_pop_settings))

        self.output = ""
        self.cur_passage = 0

    def next_passage(self):
        self.cur_passage += 1
        # handle events
        while self.cur_passage in self.settings["events"][0]:
            loc = self.settings["events"][0].index(self.cur_passage)
            self.settings[self.settings["events"][1][loc]] = self.settings["events"][2][
                loc
            ]
            del self.settings["events"][0][loc]
            del self.settings["events"][1][loc]
            del self.settings["events"][2][loc]

        for gen in range(self.settings["n_gen_per_transfer"] - 1):
            for pop in self.sims:
                pop.new_generation()

        # passage: change population size, handle migration
        previous = []
        for pop, i in zip(self.sims, list(range(self.n_pop))):
            # keep previous generation for migration
            if self.settings["migration"] is not None:
                previous.append(pop.current_population)

            # sample
            self.output += pop.current_population.to_fasta(
                n_seq=self.settings["sampling"][i],
                description="-pop{}-gen{} ".format(i, self.cur_passage),
            )

            # change parameters for transfer
            if "transfer_amount" in self.settings:
                pop.settings.max_pop = self.settings["transfer_amount"]
            else:
                pop.settings.max_pop = (
                    self.settings["transfer_prop"] * pop.current_population.n_seq
                )
            pop.new_generation()
            try:
                pop.settings.max_pop = self.settings["max_pop"][i]
            except IndexError:
                pop.settings.max_pop = self.settings["max_pop"]

        # handle migration
        if self.settings["migration"] is not None:
            for origin, target in permutations(range(self.n_pop), 2):
                n_migrate = int(
                    self.settings["migration"][target][origin] * previous[origin].n_seq
                )  # number of sequences that migrate
                sample = previous[origin].get_sample(n_migrate)  # sampled sequences
                for seq, i in zip(sample, list(range(n_migrate))):
                    changed = previous[origin].get_seq(
                        seq
                    )  # get the changes is this sequence
                    self.sims[target].current_population.add_sequence(changed)

    def all_passages(self):
        for passage in tqdm(list(range(self.settings["n_transfer"]))):
            self.next_passage()


def run(scenario_settings, organism_settings):
    passaging_run = Passaging(organism_settings, scenario_settings)
    passaging_run.all_passages()
    fasta = passaging_run.output
    if len(fasta) < 6000:
        print(passaging_run.sims)
    return fasta


def main():
    """Entry point for virus passaging."""
    # parse command line arguments
    parser = argparse.ArgumentParser(
        description="simulation of a virus passaging experiment"
    )
    parser.add_argument(
        "-n",
        "--transfer-number",
        type=int,
        default=10,
        help="number of transfers target run the simulation for",
    )
    parser.add_argument(
        "-g",
        "--generations",
        type=int,
        default=2,
        help="number of generations per transfer, defaults target 2",
    )
    parser.add_argument(
        "--organism",
        nargs=1,
        default="HIV",
        help="organism simulation settings target use",
    )
    parser.add_argument(
        "--init",
        nargs="+",
        default=[1],
        type=int,
        help="list of initial population sizes for each population, \
                              or single value (use same initial pop size for all populations),\
                              defaults target 1 per population",
    )
    parser.add_argument(
        "--max-population",
        nargs="+",
        default=[10000],
        type=int,
        help="list of maximum population sizes for each population, \
                              or single value (use same max pop size for all populations),\
                              defaults target 10000 per population",
    )
    parser.add_argument(
        "--transfer",
        nargs="+",
        default=1,
        type=float,
        help="list of transfer proportions for each population, \
                              or single value (use same transfer proportion for all populations),\
                              defaults target 0.01 per population",
    )
    parser.add_argument(
        "--migration",
        default="0",
        type=str,
        help="migration rate matrix, \
                              or single value (use migration rate between all populations),\
                              defaults target no migration",
    )
    parser.add_argument(
        "--sampling",
        nargs="+",
        default=[0],
        type=int,
        help="list of sampling sizes for each population, \
                              or single value (use same sampling size for all populations),\
                              defaults target no sampling",
    )
    parser.add_argument(
        "--events",
        default="",
        type=str,
        help='list of events (changes in parameters) during simulation, \
                            format: "[passage number] [parameter target change] [new values] : \
                            [passage number] [parameters target change] [new values] : ..."',
    )

    parser.add_argument("--output", type=str, help="output file name")

    args = parser.parse_args()
    settings = {}

    settings["n_transfer"] = args.transfer_number

    settings["pop_size"] = args.init

    settings["max_pop"] = args.max_population

    settings["n_gen_per_transfer"] = args.generations

    settings["transfer_prop"] = args.transfer

    if "[" not in args.migration:
        settings["migration"] = np.ones([len(args.init)] * 2) * float(args.migration[0])
    else:
        settings["migration"] = ast.literal_eval(args.migration)

    settings["sampling"] = args.sampling

    events = [[], [], []]
    for i in args.events.split(":"):
        if len(i) > 0:
            fields = i.split()
            events[0].append(int(fields[0]))
            events[1].append(fields[1])
            if "." in i:
                events[2].append([float(i) for i in fields[2:]])
            else:
                events[2].append([int(i) for i in fields[2:]])
    settings["events"] = events

    # run scenario
    simulation_settings = SimulationSettings.from_preset(args.organism)
    out = run(settings, simulation_settings)
    if args.output:
        with open(args.output, "w") as fd:
            fd.write(out)
    else:
        print(out)


if __name__ == "__main__":
    main()
