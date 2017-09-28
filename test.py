#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 16:24:41 2017

@author: eva
"""

import recreate_dataset as rc

sample_sizes = '33 37 25 25 49 24 39 30 23 62 42 42 50 57 36 46 45 39 25 18 11 14 17 11 31 20 26 20 18 27 21 30 47 53 40 45 27 28 35 19 13 36 26 17 25 34 22 18 23 19 23 21 18 20 39 20 22 23 26 54 42 46 32 14 19 27 24 32 21 31 29 19 40 38 36 29 32 35 27 16 28 20 67 29 42 43 16 21 15 30 26 15 25 17 10 21 18 17'.split(' ')
sample_sizes = [int(i) for i in sample_sizes] 
    
nr_mutations = '17 29 19 20 36 10 17 14 15 44 12 23 40 20 14 15 35 38 17 7 6 7 8 4 7 24 10 9 8 24 4 15 38 31 28 28 18 11 30 10 4 11 6 9 14 7 8 12 14 9 5 10 11 15 11 12 5 9 11 27 15 23 44 10 20 12 11 12 13 25 13 1 14 23 29 13 13 11 13 12 9 5 31 13 21 33 11 13 10 40 14 6 22 31 5 16 10 17'.split(' ')
nr_mutations = [int(i) for i in nr_mutations]

apobec = '1 1 0 0 1 0 0 0 0 1 0 0 1 1 0 0 1 1 1 0 0 0 0 0 0 1 0 0 1 1 0 0 0 0 1 1 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0 1 0 1 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0'.split(' ')
apobec = [int(i) for i in apobec]
    
model = 'lognormal'
parameters = {'fl':0.14,'mu':-0.22, 'sigma':0.19}
apobec_rate = 25.99
    
action = 'shared_stats'
    
simul, s, t = rc.recreate_dataset(sample_sizes[:5],nr_mutations,apobec,model,parameters,
                                   apobec_rate,action=action)

simul.ga_increase = 20
for i in range(20):
    simul.new_generation()
simul_stats = simul.current_gen.stats()