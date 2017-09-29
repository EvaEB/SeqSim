# -*- coding: utf-8 -*-
"""
Created on Mon Aug 21 11:18:40 2017

@author: evabo
"""
import numpy as np
import simulation as sim
import sys
#sys.path.append('C:\\users\\evabo\\polybox\\PhD\\code\\sequence_vis')
sys.path.append('/home/eva/polybox/PhD/code/sequence_tools')

import highlighter

class transmission_tree():
    def __init__(self,tree):
        self.tree = tree
        self.patients = tree.keys()
        self.n_gen = len(tree[tree.keys()[0]])
        self.sims = {i: None for i in self.patients}
        
        
    def __str__(self):
        string = ''
        for i in tree:
            for j in i:
                string+= j+'\t'+tree[j]
            string+='\n'
        return string
        
    def get_fasta(self,n_seq=5):
        fasta = ''
        for patient in self.patients:
            fasta += self.sims[patient].current_gen.to_fasta(
                        n_seq=n_seq, description='-{}'.format(patient))
        return fasta
        
    
    def simulate_tree(self):
        for gen in range(self.n_gen):
            for patient in self.patients:
                if self.tree[patient][gen] == '>':
                    self.sims[patient] = sim.Simulation(simulation_settings='HIV',
                                                        name=patient,seq_len=2600)
                    self.sims[patient].new_generation()
                if self.tree[patient][gen].isalpha():
                    counter=0
                    while True:
                        counter+=1
                        self.sims[patient] = self.sims[self.tree[patient][gen]].copy(name=patient,n_seq=1)
                        if self.sims[patient].get_nr_offspring(0,return_fitness=True)[1] > 0.5:
                            print counter
                            break
                if self.tree[patient][gen] == '-':
                    self.sims[patient].new_generation()
                if self.tree[patient][gen] == '=':
                    for i in range(10):
                        self.sims[patient].new_generation()
                

if __name__ == '__main__':
    #O: not present
    #-: singe generation
    #\: transmission from sequence above
    #/: transmission from sequence below
    #=: 10 generations
    tree = {'A': '........B====---====--+',
            'B': '>-=====--====-------+..',
            'C': '...............B====+..',
            'D': '..............A======-+'}
    
    #tree = np.array([list(i) for i in tree])
    t = transmission_tree(tree)
    print t 
    t.simulate_tree()
    fasta = t.get_fasta()
    seqs,labels = highlighter.fasta_to_seq(fasta)
    consensus = highlighter.consensus(seqs)
    highlighter.highlighter_plot(seqs,consensus,labels)
    
    with open('/home/eva/polybox/PhD/code/sequence_tools/test.fasta','w') as f:
        f.write(t.get_fasta())