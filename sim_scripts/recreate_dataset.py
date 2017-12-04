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



def recreate_dataset(sample_sizes,nr_mutations,apobec,model,parameters,
                     apobec_rate,action='print'):
    simul = sim.Simulation(simulation_settings='HIV',model=model,parameters=parameters)
    previous = 0
    if action == 'shared_stats':
        stats_here = stats()
    elif action == 'n_generations':
        n_gen = []
    for e, sample_size in enumerate(sample_sizes):
        this_patient = simul.copy(name=e)

        for i in range(50):
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
                        print i
                        n_gen.append(i)
                    break
                else:
                    if action == 'print':
                        print '#{} generations'.format(i-1)
                        previous_gen.print_sample(old_sample)
                    elif action == 'shared_stats':
                        stats_here.add(previous_gen,old_sample,e)
                    elif action == 'n_generations':
                        print i
                        n_gen.append(i)
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
    sample_sizes = [58] #'33 37 25 25 49 24 39 30 23 62 42 42 50 57 36 46 45 39 25 18 11 14 17 11 31 20 26 20 18 27 21 30 47 53 40 45 27 28 35 19 13 36 26 17 25 34 22 18 23 19 23 21 18 20 39 20 22 23 26 54 42 46 32 14 19 27 24 32 21 31 29 19 40 38 36 29 32 35 27 16 28 20 67 29 42 43 16 21 15 30 26 15 25 17 10 21 18 17'.split(' ')
    sample_sizes = [int(i) for i in sample_sizes]

    nr_mutations = [20] #'17 29 19 20 36 10 17 14 15 44 12 23 40 20 14 15 35 38 17 7 6 7 8 4 7 24 10 9 8 24 4 15 38 31 28 28 18 11 30 10 4 11 6 9 14 7 8 12 14 9 5 10 11 15 11 12 5 9 11 27 15 23 44 10 20 12 11 12 13 25 13 1 14 23 29 13 13 11 13 12 9 5 31 13 21 33 11 13 10 40 14 6 22 31 5 16 10 17'.split(' ')
    nr_mutations = [int(i) for i in nr_mutations]

    apobec = [1]#'1 1 0 0 1 0 0 0 0 1 0 0 1 1 0 0 1 1 1 0 0 0 0 0 0 1 0 0 1 1 0 0 0 0 1 1 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0 1 0 1 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0'.split(' ')
    apobec = [int(i) for i in apobec]

    model = 'exponential'
    parameters = {'fb':0.2, 'lb':0.6, 'fn':0.8}
    #apobec_rate = np.exp(1.830)
    #model = 'lognormal'
    #parameters = {'fl':0.045,'mu':-0.248, 'sigma':0.149}
    #model = 'neutral'
    #parameters = {}
    apobec_rate = 10.621

    action = 'n_generations'

    for rep in range(1):
        print rep
        times_since_infection  = recreate_dataset(sample_sizes,nr_mutations,apobec,model,parameters,
                                              apobec_rate,action=action)

        print '12 (8-16) \t 19 (9, 28)'
        #pickle.dump(stats_final,open('/home/eva/polybox/PhD/FED/SMC/best_fit_sims/lognormal_new/stats_{}.npy'.format(rep),'w'))
