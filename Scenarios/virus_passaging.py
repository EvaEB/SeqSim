"""
Created on Tue Jul 11 15:13:13 2017

@author: eva

virus passaging as the experiment in the liquid handling robot
"""
from __future__ import print_function
import os
import sys
from copy import deepcopy as copy

dir_path = os.path.dirname(os.path.realpath(__file__))
dir_path_up = os.sep.join(dir_path.split(os.sep)[:-1])
sys.path.append(dir_path_up+os.sep+'SeqSimEvo')
import SeqSimEvo as sim

from tqdm import tqdm
import yaml


class passaging():
    def __init__(self, sim_settings, passaging_settings):
        self.settings = copy(passaging_settings)
        main_sim = sim.Simulation(copy(sim_settings), n_seq_init = max(self.settings['pop_size']))
        self.sims = []
        self.n_pop = len(self.settings['pop_size'])

        names = iter(range(self.n_pop))
        for i in range(self.n_pop):
            self.sims.append(main_sim.copy(next(names),n_seq=self.settings['pop_size'][i]))
            try:
                self.sims[i].settings['max_pop'] = self.settings['max_pop'][i]
            except TypeError:
                self.sims[i].settings['max_pop'] = self.settings['max_pop']
        self.output = main_sim.current_gen.to_fasta(n_seq=1,description='consensus')
        self.cur_passage = 0

    def next_passage(self):
        self.cur_passage+=1
        #handle events
        while self.cur_passage in self.settings['events'][0]:
            loc = self.settings['events'][0].index(self.cur_passage)
            self.settings[self.settings['events'][1][loc]] = self.settings['events'][2][loc]
            del self.settings['events'][0][loc]
            del self.settings['events'][1][loc]
            del self.settings['events'][2][loc]

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

            self.output+= pop.current_gen.to_fasta(n_seq=self.settings['sampling'][i],
                                                   description=' - pop {} - transfer {} '.format(i,self.cur_passage))

            #change parameters for transfer
            if 'transfer_amount' in self.settings:
                pop.settings['max_pop'] = self.settings['transfer_amount']
            else:
                pop.settings['max_pop'] = self.settings['transfer_prop']*pop.current_gen.n_seq
            pop.new_generation()
            try:
                pop.settings['max_pop'] = self.settings['max_pop'][i]
            except TypeError:
                pop.settings['max_pop'] = self.settings['max_pop']

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

    def all_passages(self):
        for passage in tqdm(range(self.settings['n_transfer'])):
            self.next_passage()


def run(scenario_settings,organism_settings):
    passaging_run = passaging(organism_settings,scenario_settings)
    passaging_run.all_passages()
    fasta = passaging_run.output
    if len(fasta) < 6000:
        print(passaging_run.sims)
    return fasta



if __name__ == '__main__':
    import argparse
    import numpy as np
    import ast
    #parse command line arguments
    parser = argparse.ArgumentParser(description='simulation of a virus passaging experiment')
    parser.add_argument('-n', type=int,default=10,
                        help='number of transfers to run the simulation for')
    parser.add_argument('-gen',type=int,default=2,
                        help='number of generations per transfer, defaults to 2')
    parser.add_argument('-o', nargs=1,default = 'HIV',
                        help='organism simulation settings to use')
    parser.add_argument('-init',nargs='+',default=[1],type=int,
                        help='list of initial population sizes for each population, \
                              or single value (use same initial pop size for all populations),\
                              defaults to 1 per population')
    parser.add_argument('-pop',nargs='+',default=[10000],type=int,
                        help='list of maximum population sizes for each population, \
                              or single value (use same max pop size for all populations),\
                              defaults to 10000 per population')
    parser.add_argument('-transfer',nargs='+',default=1,type=float,
                        help='list of transfer proportions for each population, \
                              or single value (use same transfer proportion for all populations),\
                              defaults to 0.01 per population')
    parser.add_argument('-mig',default='[[0,0],[0,0]]',type=str,
                        help='migration rate matrix, \
                              or single value (use migration rate between all populations),\
                              defaults to no migration')
    parser.add_argument('-sampling',nargs='+',default=[0],type=int,
                        help='list of sampling sizes for each population, \
                              or single value (use same sampling size for all populations),\
                              defaults to no sampling')
    parser.add_argument('-events',default='',type=str,
                        help='list of events (changes in parameters) during simulation, \
                            format: "[passage number] [parameter to change] [new values] : \
                            [passage number] [parameters to change] [new values] : ..."')

    args = parser.parse_args()
    settings = {}

    settings['n_transfer'] = args.n

    settings['pop_size'] = args.init

    settings['max_pop'] = args.pop

    settings['n_gen_per_transfer'] = args.gen

    settings['transfer_prop'] = args.transfer

    if '[' not in args.mig:
        settings['migration'] = np.ones([len(args.init)]*2)*float(args.mig[0])
    else:
        settings['migration'] = ast.literal_eval(args.mig)


    settings['sampling'] = args.sampling

    events = [[],[],[]]
    for i in args.events.split(':'):
        if len(i) > 0:
            fields = i.split()
            events[0].append(int(fields[0]))
            events[1].append(fields[1])
            if '.' in i:
                events[2].append([float(i) for i in fields[2:]])
            else:
                events[2].append([int(i) for i in fields[2:]])
    settings['events'] = events


    #run scenario
    print(run(settings, args.o))
