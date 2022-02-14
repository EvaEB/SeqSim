"""Simulate robot experiment as executed by peyolution."""

import os
import argparse

import pandas as pd

from peyolution.transmissions import transmissions

from SeqSimEvo import Sequence, Simulation, SimulationSettings, Population


def main(args):
    """Main function."""
    simulation_settings = SimulationSettings.from_preset(
        args.preset, dilution=args.dilution
    )
    sequence = Sequence.generate_sequence(simulation_settings.seq_len)
    compartments = [
        Simulation(sequence, simulation_settings, name=compartment)
        for compartment in range(args.n_compartments)
    ]
    plan = pd.read_csv(args.passage_plan)
    for generation in range(args.generations):
        print(f"{generation=}")
        events = plan[plan.generation == generation]

        if len(events) <= 0:
            for compartment in compartments:
                compartment.new_generation()
            continue

        transfer = None
        for _, event in events.iterrows():
            if event.event == "transmission":
                transfer = transmissions[event.value]
                break

        if transfer is None:
            continue

        splits = [[] for _ in range(args.n_compartments)]
        for step in transfer:
            for origin, volume in zip(step.origin, step.volumes):
                splits[origin // 2] += [volume]

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

        populations = [[] for _ in range(args.n_compartments)]
        for idx, step in enumerate(transfer):
            for origin, target in zip(step.origin, step.target):
                populations[target // 2] += [split_populations[origin // 2][idx]]

        merged_populations = [Population.merge(*pop) for pop in populations]

        for compartment, population in zip(compartments, merged_populations):
            compartment.current_population = population
            compartment.gen += 1

    output = "".join(
        [
            compartment.current_population.to_fasta(n_seq=1, description=f"pop-{idx}")
            for idx, compartment in enumerate(compartments)
        ]
    )
    if args.output:
        if not os.path.exists(os.path.dirname(args.output)):
            os.makedirs(os.path.dirname(args.output))
        with open(args.output, "w", encoding="utf8") as file_descriptor:
            file_descriptor.write(output)
    else:
        print(output)


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
        "--dilution",
        type=float,
        help="Dilution at every generation.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        help="Output file name. If not set, return to stdout.",
    )

    _args = parser.parse_args()
    main(_args)
