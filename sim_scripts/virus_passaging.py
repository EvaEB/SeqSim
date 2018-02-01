#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 11 15:13:13 2017

@author: eva

virus passaging as the experiment in the liquid handling robot
"""
#import optparse
import seq_sim as sim
#from tqdm import tqdm
import yaml

class passaging():
    def __init__(self, sim_settings, passaging_settings):
        self.settings = passaging_settings

        main_sim = sim.Simulation(sim_settings, n_seq_init = max(self.settings['pop_size']))
        self.sims = []
        self.n_pop = len(self.settings['pop_size'])

        names = iter(range(self.n_pop))
        for i in range(self.n_pop):
            self.sims.append(main_sim.copy(names.next(),n_seq=self.settings['pop_size'][i]))
            self.sims[i].settings['max_pop'] = self.settings['max_pop'][i]
        self.output = ''
        self.cur_passage = 0

    def next_passage(self):
        self.cur_passage+=1
        #handle events
        while self.cur_passage in self.settings['events'][0]:
            loc = self.settings['events'][0].index(self.cur_passage)
            self.settings[self.settings['events'][1][loc]] = self.settings['events'][2][loc]
            del self.settings['events'][0][loc]

        for gen in range(self.settings['n_gen_per_transfer']-1):
            for pop in self.sims:
                pop.new_generation()


        #passage: change population size, handle migration
        previous = []
        for pop,i in zip(self.sims, range(self.n_pop)):
            #keep previous generation for migration
            if self.settings['migration'] is not None:
                previous.append(pop.current_gen)

            #sample
            self.output+= pop.current_gen.to_fasta(n_seq=self.settings['sampling'][i])

            #change parameters for transfer
            pop.settings['max_pop'] = self.settings['transfer_prop']*pop.current_gen.n_seq
            pop.new_generation()
            pop.settings['max_pop'] = self.settings['max_pop'][i]

        #handle migration
        if self.settings['migration'] is not None:
            changed_index = []
            for to in range(self.n_pop):
                for fro in range(self.n_pop):
                    if to != fro:
                        n_migrate = int(self.settings['migration'][to][fro]*previous[fro].n_seq) #number of sequences that migrate
                        sample = previous[fro].get_sample(n_migrate) #sampled sequences
                        for seq,i in zip(sample,range(n_migrate)):
                            changed = previous[fro].get_seq(seq) #get the changes is this sequence
                            self.sims[to].current_gen.add_sequence(changed)

if __name__ == '__main__':
    with open('/home/eva/code/SeqSim/seq_sim/simulation_settings/small') as sim_settings:
        with open('/home/eva/code/SeqSim/sim_scripts/settings_files/RecreateDataset_template') as pas_settings:
            passaging = passaging(yaml.safe_load(sim_settings),yaml.safe_load(pas_settings))

    for i in range(3):
        passaging.next_passage()
