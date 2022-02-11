"""
recreate_dataset.py: a SeqSimEvo scenario

Recreate a dataset of sequence samples from different patients. The amount of
unique mutations per sequence sample is matched as closely as possible. A version
of this code was used for https://doi.org/10.1093/ve/vey029
"""

from collections import Counter
import numpy as np
import pickle
from SeqSimEvo import simulation as sim


class stats(object):
    def __init__(self):
        self.mutation_counts = {}
        self.n_patients = -1

    def add(self, population, sample, name):
        self.n_patients += 1
        all_muts = []
        for i in sample:
            changes = population.get_seq(i)
            if changes is not None:
                for mut in changes:
                    all_muts.append(str(mut))
        for mut in set(all_muts):
            try:
                self.mutation_counts[mut].append(self.n_patients)
            except KeyError:
                self.mutation_counts[mut] = [self.n_patients]

    def get_stats(self, sim):
        values = list(self.mutation_counts.values())
        shared_index_temp = np.array(
            [[patient, len(mutlist)] for mutlist in values for patient in mutlist]
        )
        shared_index = [
            np.mean(shared_index_temp[shared_index_temp[:, 0] == i + 1, 1] - 1)
            for i in range(self.n_patients)
        ]

        shared_histogram = dict(Counter([len(i) for i in values]))

        fitnesses = {}
        shared_fitness = []
        for mut_id in list(self.mutation_counts.keys()):
            mut = mut_id.strip("[]").split()
            shared_fitness.append(
                [
                    len(self.mutation_counts[mut_id]),
                    sim.get_fitness_effect(int(mut[0]), int(mut[1])),
                ]
            )
        shared_fitness = np.array(shared_fitness)
        multiplicities = set(shared_fitness[:, 0])
        for i in multiplicities:
            fitnesses[i] = shared_fitness[shared_fitness[:, 0] == i, 1]

        return shared_index, shared_histogram, fitnesses


def recreate_dataset(
    sample_sizes,
    nr_mutations,
    apobec,
    model=None,
    parameters=None,
    apobec_rate=1,
    action="print",
    simulation_settings="HIV",
):
    if model is not None:
        simul = sim.Simulation(simulation_settings, model=model, parameters=parameters)
    else:
        simul = sim.Simulation(simulation_settings)
    previous = 0
    if action == "shared_stats":
        stats_here = stats()
    elif action == "n_generations":
        n_gen = []
    elif action == "fasta":
        fasta = simul.current_gen.to_fasta(n_seq=1, description=" - ancestor")
    for e, sample_size in enumerate(sample_sizes):
        print(e)
        this_patient = simul.copy(name=e)
        for i in range(150):
            this_patient.new_generation()
            if apobec[e] == 1:
                this_patient.ga_increase = apobec_rate
            else:
                this_patient.ga_increase = 1
            sample = this_patient.current_gen.get_sample(sample_size)
            changes = []
            for seq_id in sample:
                try:
                    changes += [
                        str(mut) for mut in this_patient.current_gen.get_seq(seq_id)
                    ]
                except TypeError:
                    pass
            n_changed = len(set(changes))
            if n_changed > nr_mutations[e]:
                if abs(nr_mutations[e] - n_changed) < abs(nr_mutations[e] - previous):
                    if action == "print":
                        print("#{} generations".format(i))
                        this_patient.current_gen.print_sample(sample)
                    elif action == "shared_stats":
                        stats_here.add(this_patient.current_gen, sample, e)
                    elif action == "n_generations":
                        # print i
                        poisson = this_patient.current_gen.Hamming_distance(
                            this_patient.settings, sample, action="print"
                        )
                        n_gen.append([i, poisson])
                    elif action == "fasta":
                        fasta += this_patient.current_gen.to_fasta(
                            seq_ids=sample,
                            description=" - patient {}".format(
                                this_patient.settings["name"]
                            ),
                        )
                    break
                else:
                    if action == "print":
                        print("#{} generations".format(i - 1))
                        previous_gen.print_sample(old_sample)
                    elif action == "shared_stats":
                        stats_here.add(previous_gen, old_sample, e)
                    elif action == "n_generations":
                        # print i
                        poisson = previous_gen.Hamming_distance(
                            this_patient.settings, sample, action="print"
                        )
                        n_gen.append([i, poisson])
                    elif action == "fasta":
                        fasta += previous_gen.to_fasta(
                            seq_ids=sample,
                            description=" - patient {}".format(
                                this_patient.settings["name"]
                            ),
                        )
                    break
            else:
                previous = n_changed
                old_sample = sample
                previous_gen = this_patient.current_gen.copy()
    if action == "shared_stats":
        return this_patient, stats_here, stats_here.get_stats(simul)
    elif action == "n_generations":
        return n_gen
    elif action == "fasta":
        return fasta


def main():
    import argparse

    # parse command line arguments
    parser = argparse.ArgumentParser(
        description="recreate a dataset of sequence samples with a given amount of mutations and sequences per sample"
    )
    parser.add_argument(
        "-o", nargs=1, default=["HIV"], help="organism simulation settings to use"
    )
    parser.add_argument(
        "-p",
        nargs=1,
        default=None,
        help="path to file with samplesizes, number of mutations and apobec status",
    )
    parser.add_argument(
        "-samples",
        nargs="+",
        default=None,
        type=int,
        help="number of samples for each individual to recreate",
    )
    parser.add_argument(
        "-muts",
        nargs="+",
        default=None,
        type=int,
        help="number of unique mutations in each individual to recreate",
    )
    parser.add_argument(
        "-apobec",
        nargs="+",
        default=None,
        type=int,
        help="apobec status for all of the individuals to recreate",
    )
    parser.add_argument(
        "-action",
        default="print",
        help="action: print: print summary of changes\
                             fasta: output fasta,\
                             shared_stats: some statistics on sharing,\
                             n_generations: the number of generations simulated",
    )
    parser.add_argument(
        "-apobec_rate",
        default=1,
        type=int,
        help="change in G-to-A mutation rate due to APOBEC",
    )
    args = parser.parse_args()

    settings = args.o[0]

    action = args.action
    apobec_rate = args.apobec_rate

    if args.p is not None:
        patient_file = args.p[0]
        with open(patient_file) as f:
            patients = f.readlines()

        nr_mutations = patients[1].strip().split(" ")
        nr_mutations = [int(i) for i in nr_mutations]

        sample_sizes = patients[3].strip().split(" ")
        sample_sizes = [int(i) for i in sample_sizes]

        apobec = patients[5].strip().split(" ")
        apobec = [int(i) for i in apobec]
    else:
        if args.samples is not None:
            sample_sizes = args.samples
        else:
            sample_sizes = [10]
        if args.muts is not None:
            nr_mutations = args.muts
        else:
            nr_mutations = [10]

        if args.apobec is not None:
            apobec = args.apobec
        else:
            apobec = [0]

    simulation = recreate_dataset(
        sample_sizes,
        nr_mutations,
        apobec,
        apobec_rate=apobec_rate,
        action=action,
        simulation_settings=settings,
    )
    if action == "shared_stats":  # todo: prettier print
        for i in simulation:
            print(i)

    if (action == "n_generations") or (action == "fasta"):
        print(simulation)


if __name__ == "__main__":
    main()
