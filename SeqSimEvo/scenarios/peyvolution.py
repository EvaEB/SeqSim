"""Simulate robot experiment as executed by pyolution."""

import argparse
from itertools import permutations

import pandas as pd
import numpy as np

from peyolution.transmissions import transmissions

from SeqSimEvo import Sequence, Simulation, SimulationSettings, Population


def main(args):
    """Main function."""
    simulation_settings = SimulationSettings.from_preset(args.preset)
    sequence = Sequence.generate_sequence(simulation_settings.seq_len)
    compartments = [
        Simulation(sequence, simulation_settings, name=compartment)
        for compartment in args.n_compartments
    ]
    plan = pd.read_csv(args.passage_plan)
    for generation in range(args.generations):
        events = plan[plan.generation == generation]

        if event <= 0:
            for compartment in compartments:
                compartment.new_generation()
            continue

        transfer = None
        for _, event in events:
            if event.event == "transmisson":
                transfer = transmissions[event.value]
                break
        if transfer is None:
            continue

        splits = [[]] * args.n_compartments
        for step in transfer:
            for origin, volume in zip(step.origin, step.volumes):
                splits[origin] += [volume]

        compartment_weights = [
            compartment.generate_offspring() for compartment in compartments
        ]

        compartment_offspring = [
            compartment.subsample_population(weights, sum(split), sort=False)
            for compartment, weights, split in zip(
                compartments, compartment_weights, splits
            )
        ]

        compartment_populations = [
            compartment.population_from_offspring(offspring)
            for compartment, offspring in zip(compartments, compartment_offspring)
        ]

        split_populations = [
            population.split(split)
            for population, split in zip(compartment_populations, splits)
        ]

        populations = [[]] * args.n_compartments
        for idx, step in enumerate(transfer):
            for origin, target in zip(step.origin, step.target):
                populations[target] += split_populations[origin][idx]

        merged_populations = [Population.merge(pop) for pop in populations]

        for compartment, population in zip(compartments, merged_populations):
            compartment.current_population = population
            compartment.gen += 1


def entry():
    """SeqSimEvo Pyolve entry point."""
    parser = argparse.ArgumentParser(description="Simulate robot experiment.")

    parser.add_argument(
        "-n",
        "--n-compartments",
        type=int,
        default=4,
        help="Number of compartments.",
    )
    parser.add_argument(
        "-g",
        "--generations",
        type=int,
        default=200,
        help="Number of transfers.",
    )
    parser.add_argument(
        "--preset",
        type=str,
        default="phix174",
        help="Name of preset used for the simulation settings.",
    )
    parser.add_argument(
        "--init-size",
        default=1e5,
        type=int,
        help="Initial population size.",
    )
    parser.add_argument(
        "--max-population",
        default=1e9,
        type=int,
        help="Maximum population size to be reached.",
    )
    parser.add_argument(
        "--passage-plan",
        type=str,
        help="Path to passage plan that includes all events during passaging.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        help="Output file name. If not set, return to stdout.",
    )

    _args = parser.parse_args()
    main(_args)