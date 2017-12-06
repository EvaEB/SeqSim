#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 11 15:13:13 2017

@author: eva
"""
import optparse
import seq_sim as sim
from tqdm import tqdm


# def passaging(settings,n_gen=50,initial_size=1,max_size=1e5,transfer_size=0.005,
#               transfer_time=2):
#
#     phage_sim = sim.Simulation(simulation_settings=settings,
#                                n_seq_init=initial_size,
#                                max_pop = max_size)
#
#
#
#     for gen in range(1,n_gen+1):
#         if gen % transfer_time == 0:
#             transferred = int(transfer_size*phage_sim.current_gen.n_seq)
#             if transferred > 1:
#                 phage_sim.settings['max_pop'] = transferred
#             else:
#                 phage_sim.settings['max_pop'] = 2
#             phage_sim.new_generation()
#             phage_sim.settings['max_pop'] = max_size
#         else:
#             phage_sim.new_generation()
#
#     return phage_sim

class passaging(sim.Simulation):
    def __init__(self,simulation_settings,initial_size,max_size,
                            transfer_prop,transfer_time, **kwargs):
        sim.Simulation.__init__(self,simulation_settings,n_seq_init=initial_size,**kwargs)
        self.initial_size = initial_size
        self.max_size = max_size
        self.transfer_prop = transfer_prop
        self.transfer_time = transfer_time


    def passage(self,n_passage,progressbar=False):
        if progressbar:
            bar = tqdm(range(n_passage))

        for i in range(n_passage):
            for j in range(self.transfer_time-1):
                self.new_generation()
                if progressbar:
                    bar.update()
            self.settings['max_pop'] = max(2,self.transfer_prop*self.current_gen.n_seq)
            self.new_generation()
            self.settings['max_pop'] = self.max_size
            
        if progressbar:
            bar.close()


if __name__ ==  '__main__':
    parser = optparse.OptionParser()

    parser.add_option("-g", "--n_gen", dest="n_gen",
                      help="number of generations to run the simulation for")

    parser.add_option("-p", "--passages", dest="n_pass",
                      help="number of passages to run the simulation for")

    parser.add_option("-i", "--initial_size", dest="initial_size",
                      help="initial size of the virus population (default 1)")

    parser.add_option("-m", "--max_size", dest="max_size",
                      help="maximum size of the virus population")

    parser.add_option("-t", "--transfer_time", dest="transfer_time",
                      help="number of generations before each transfer")

    parser.add_option("-r", "--transfer_prop", dest="transfer_prop",
                      help="proportion of viral population that is transferred")

    parser.add_option("-s", "--settings", dest="settings",
                      help="settings file for the simulation")

    (options, args) = parser.parse_args()

    if options.transfer_time is None:
        options.transfer_time = 2

    if options.n_gen is None and options.n_pass is None:
        parser.error('either number of generations (-g) or number of passages (-p) is required')
    elif options.n_gen is not None and options.n_pass is not None:
        parser.error('both number of generations (-g) or number of passages (-p) provided, please only provide one of the two')
    elif options.n_pass is not None:
        options.n_gen = options.n_pass * options.transfer_time

    if options.settings is None:
        options.settings = 'phix174'
    if options.max_size is None:
        options.max_size = 1e5

    if options.transfer_prop is None:
        options.transfer_prop = 0.005

    if options.initial_size is None:
        options.initial_size = 1e5

    phage_sim = passaging(options.settings,options.initial_size,options.max_size,
                          float(options.transfer_prop),options.transfer_time)

    phage_sim.passage(100,progressbar=True)

    print phage_sim
