#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 24 17:05:10 2017

@author: eva
"""

import simulation as simulation
import random

class patient(object):
    def __init__(self,sim,donor='>',name=0,gen=0):
        self.name = name
        self.sim = sim
        self.active=True
        self.donor = donor
        self.creation_time = gen
        self.death_time = 0
    def __str__(self):
        return '{}>{} ({}-{})'.format(self.donor,self.name,self.creation_time,
                                      self.death_time)

class transmission_simulation(object):
    def __init__(self, sim, n_gen = 50, p_transmission = 0.2,deathrate=0.2):
        self.patients = [0]
        self.patients[0] = patient(sim)
        self.n_gen = n_gen
        self.p_transmission = p_transmission
        self.sim = sim
        self.death_rate = 0.2
        
    def simulate_transmission(self):
        for gen in range(self.n_gen):
            for pat in self.patients:
                if pat.active:
                    if random.random() < self.p_transmission:
                        name = len(self.patients)
                        self.patients.append(patient(self.sim.copy(name=name),
                                                     donor=pat.name,
                                                     name=len(self.patients),
                                                     gen = gen))
                    if random.random() < self.death_rate:
                        pat.active = False
                        pat.death_time = gen
                    else:
                        pat.sim.new_generation()

    def __str__(self):
        string = ''
        for p in self.patients:
            string += '.'*p.creation_time + str(p.donor) 
            string += '-'*(p.death_time-p.creation_time)
            string += 'x'
            string += '.'*(self.n_gen-p.death_time)
            string += '\n'
        return string
                
        
if __name__ == '__main__':
    base_sim = simulation.Simulation()
    t = transmission_simulation(base_sim)
    t.simulate_transmission()
    print t
    for i in t.patients:
        print i.sim.current_gen.to_fasta(n_seq=3,name='patient'+i.name)