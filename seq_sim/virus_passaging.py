#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 11 15:13:13 2017

@author: eva
"""
import optparse
import simulation as sim


def passaging():
    n_gen = 50

    initial_size = 1e5
    max_size=1e9

    transfer_size = .005
    transfer_time = 2
    
    phage_sim = sim.Simulation(simulation_settings='phix174',
                               n_seq_init=initial_size,
                               max_pop = max_size)

    for gen in range(1,n_gen+1):
        if gen % transfer_time == 0:
            transferred = int(transfer_size * phage_sim.current_gen.n_seq)
            phage_sim.settings['max_pop'] = transferred
            phage_sim.new_generation()
            phage_sim.settings['max_pop'] = max_size
        else:
            phage_sim.new_generation()
    print 'seq now, ',phage_sim.current_gen.n_seq
        
#print phage_sim
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
    

    if options.n_gen is None and options.n_pass is None:
        parser.error('either number of generations (-g) or number of passages (-p) is required')
    if options.n_gen is not None and options.n_pass is not None:
        parser.error('both number of generations (-g) or number of passages (-p) provided, please only provide one of the two')
    
    
    if options.initial_size is None:
        options.initial_size = 1
        
    print options.initial_size
    
