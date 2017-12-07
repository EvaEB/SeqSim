#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 14:29:00 2017

@author: eva
"""

import seq_sim as sim
from collections import Counter
import numpy as np
import pickle

class stats(object):
    def __init__(self):
        self.mutation_counts = {}
        self.n_patients = -1

    def add(self,population,sample,name):
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

    def get_stats(self,sim):
        values = self.mutation_counts.values()
        shared_index_temp = np.array([[patient, len(mutlist)] for mutlist in values for patient in mutlist])
        shared_index = [np.mean(shared_index_temp[shared_index_temp[:,0]==i+1,1]-1) for i in range(self.n_patients)]

        shared_histogram = dict(Counter([len(i) for i in values]))

        fitnesses = {}
        shared_fitness = []
        for mut_id in self.mutation_counts.keys():
            mut = mut_id.strip('[]').split()
            shared_fitness.append([len(self.mutation_counts[mut_id]),
                                   sim.get_fitness_effect(int(mut[0]),int(mut[1]))])
        shared_fitness = np.array(shared_fitness)
        multiplicities = set(shared_fitness[:,0])
        for i in multiplicities:
            fitnesses[i] = shared_fitness[shared_fitness[:,0]==i,1]


        return shared_index,shared_histogram,fitnesses



def recreate_dataset(sample_sizes,nr_mutations,apobec,model=None,parameters=None,
                     apobec_rate=1,action='print',simulation_settings='HIV'):
    if model is not None:
        simul = sim.Simulation(simulation_settings,model=model,parameters=parameters)
    else:
        simul = sim.Simulation(simulation_settings)
    previous = 0
    if action == 'shared_stats':
        stats_here = stats()
    elif action == 'n_generations':
        n_gen = []
    for e, sample_size in enumerate(sample_sizes):
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
                    changes+= [str(mut) for mut in this_patient.current_gen.get_seq(seq_id)]
                except TypeError:
                    pass
            n_changed = len(set(changes))
            if n_changed > nr_mutations[e]:
                if abs(nr_mutations[e]-n_changed)<abs(nr_mutations[e]-previous):
                    if action == 'print':
                        print '#{} generations'.format(i)
                        this_patient.current_gen.print_sample(sample)
                    elif action == 'shared_stats':
                        stats_here.add(this_patient.current_gen,sample,e)
                    elif action == 'n_generations':
                        #print i
                        poisson = this_patient.current_gen.Hamming_distance(this_patient.settings,sample,action='Poisson_fit')
                        n_gen.append([i, poisson])
                    break
                else:
                    if action == 'print':
                        print '#{} generations'.format(i-1)
                        previous_gen.print_sample(old_sample)
                    elif action == 'shared_stats':
                        stats_here.add(previous_gen,old_sample,e)
                    elif action == 'n_generations':
                        #print i
                        poisson = previous_gen.Hamming_distance(this_patient.settings,sample,action='Poisson_fit')
                        n_gen.append([i,poisson])
                    break
            else:
                previous = n_changed
                old_sample = sample
                previous_gen = this_patient.current_gen.copy()
    if action == 'shared_stats':
        return this_patient,stats_here, stats_here.get_stats(simul)
    elif action == 'n_generations':
        return n_gen


if __name__ == '__main__':
    import sys
    settings = sys.argv[1]
    patient_file = sys.argv[2]
    action = sys.argv[3] #print, shared_stats or n_generations
    apobec_rate = sys.argv[4]

    with open(patient_file) as f:
        patients = f.readlines()

    nr_mutations = patients[1].strip().split(' ')
    nr_mutations = [int(i) for i in nr_mutations]

    sample_sizes = patients[3].strip().split(' ')
    sample_sizes = [int(i) for i in sample_sizes]

    apobec = patients[5].strip().split(' ')
    apobec = [int(i) for i in apobec]

    simulation  = recreate_dataset(sample_sizes,nr_mutations,apobec,
                                   apobec_rate=apobec_rate,
                                   action=action,
                                   simulation_settings=settings)
    if action == 'shared_stats': #todo: prettier print
        for i in simulation:
            print i
    if action == 'n_generations':
        print simulation
