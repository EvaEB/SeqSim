#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 13:42:33 2016

@author: eva
"""
import numpy as np
np.set_printoptions(threshold=np.nan, suppress=True)
from collections import Counter
import random
from copy import copy, deepcopy
import time
import yaml
import os
import inspect
import progressbar

class Seq(object):
    '''
    Create a random sequence.

    Keyword arguments:

    * seq_len: int, length of the sequence to simulate

    * base_dist: list of length 4, base distribution in order A-G-T-C
    '''
    def __init__(self, seq_len=100, base_dist=None, seq=''):
        self.translation = 'AGTC'
        if seq == '':
            self.generate_seq(seq_len, base_dist)
            self.len = seq_len
        else:
            self.sequence = self.translate_sequence(seq)
            self.len = len(seq)



    def __str__(self):
        seq = ''
        for i in self.sequence:
            seq += self.translation[i]
        return seq

    def __len__(self):
        return self.len

    def __getitem__(self, key):
        return self.sequence[key]

    def translate_sequence(self, seq):
        '''translate a sequence from bases to numbers'''
        sequence = []
        for i in seq:
            sequence.append(self.translation.index(i))
        return sequence

    def generate_seq(self, seq_len, base_dist=None):
        '''generate a sequence of seq_len bases according to the base
        distribution (order: A-G-T-C)'''
        seq = []
        if base_dist is None:
            base_dist = np.array([0.25, 0.5, 0.75, 1])
        else:
            base_dist = np.array(base_dist)
        for i in range(seq_len):
            seq.append(sum(base_dist < random.random()))
        self.sequence = seq

class Simulation(object):
    '''
    Sequence simulation object

    initialisation parameters:
        * **model**: model to use for the MFED (otions: 'neutral'(default),\
        'lognormal','exponential','spikes','beta', 'from_data','steps')
        * **parameters**: dict containing appropriate parameters for the MFED
        * **mut_rate**: mutation rate per site per generation (default = 2.16e-5)
        * **subs_matrix**:
        * **seq_len**: length of sequence to simulate (default 2600), used if \
        no sequence is provided
        * **basedist**: list of cumulutive distribution of bases used for \
        generating the sequence in order [A,G,T,C]. (default: [0.25,0.5,0.75,1])
        * **R0**: initial average amount of offspring per sequence
        * **ga_increase**: Increase in G-A mutation rate (defaul: 1 (no increase))
        * **max_pop**: maximum population size
        * **name**
        * **sequence**: a sequence object
        * **fitness_table**
        * **n_seq_init**
    '''

    def __init__(self, simulation_settings='HIV', **kwargs):

        try:
            with open(simulation_settings) as f:
                self.settings = yaml.safe_load(f)
        except IOError:
            path = os.path.dirname(os.path.abspath(inspect.stack()[0][1]))
            with open(path+'/simulation_settings/'+simulation_settings) as f:
                self.settings = yaml.safe_load(f)
        except TypeError:
            self.settings = simulation_settings

        for key, value in kwargs.iteritems():
            self.settings[key] = value



        self.settings['subs_matrix'] = np.array(self.settings['subs_matrix'])

        #adapt mutation rate for APOBEC increase in G-to-A mutations
        self.mut_rate = sum([(self.settings['subs_matrix'][1][1]*self.settings['mut_rate']),
                             ((1-self.settings['subs_matrix'][1][1])*self.settings['ga_increase']*
                              self.settings['mut_rate'])])

        if 'sequence' not in self.settings.keys():
            self.sequence = Seq(self.settings['seq_len'], self.settings['basedist'])
        else:
            self.sequence = self.settings['sequence']
            self.settings['seq_len'] = self.sequence.len

        if 'fitness_table' not in self.settings.keys():
            if self.settings['model'] == 'exponential':
                for par in ['fl','fd','fb']:
                    if par not in self.settings['parameters'].keys():
                        self.settings['parameters'][par] = 0
                if 'fn' not in self.settings['parameters'].keys():
                    self.settings['parameters']['fn'] = 1 - (self.settings['parameters']['fl']+ \
                                                             self.settings['parameters']['fd']+ \
                                                             self.settings['parameters']['fb'])
            self.fitness_table = self.__getfitness_table()
        else:
            self.fitness_table = self.settings['fitness_table']

        self.settings['max_pop'] = int(self.settings['max_pop'])

        self.gen = 0
        self.average_fitness= 1
        self.effective_pop = 1
        self.n_seq = self.settings['n_seq_init']

        if 'pop' not in self.settings.keys():
            self.current_gen = Population(self, n_seq=int(self.settings['n_seq_init']))
        else:
            self.current_gen = pop
        #self.current_gen = np.array([[0,np.nan,np.nan,0],
        #                            [1,150,3,0],
        #                            [1,1001,0,0],
        #                            [2,20,2,0]])


        self.mutations_per_seq = self.new_mutations_per_seq()


    def new_mutations_per_seq(self):
        return iter(np.random.binomial(self.sequence.len,self.settings['mut_rate'],1000))

    def __str__(self):
        string = 'sequence simulation\n'
        string += 'MED model\t'+self.settings['model'] + '\n'
        for i in self.settings['parameters']:
            string += str(i)+'\t'+str(self.settings['parameters'][i])+'\n'
        string += 'ancestor\t'+str(self.sequence)+'\n'
        string += 'number of generations\t'+str(self.gen)+'\n'
        stats = self.current_gen.stats()
        for i in stats.keys():
            string += i.replace('_', ' ')+'\t'+str(stats[i])+'\n'
        return string


    def __getfitness_table(self):
        '''creates a table with random fitness values according to model and
        parameters par'''
        seq_len = self.sequence.len
        #create an array filled with 0 (all lethal)
        fitness = np.zeros((4, seq_len))
        #make sure the fitness benefit of the initial sequence is 1
        for (i, base) in zip(range(seq_len), self.sequence.sequence):
            fitness[base, i] = 1
        to_fill = np.where(fitness == 0)

        if self.settings['model'] == 'neutral': #neutral model
            fitness = np.ones((4, seq_len))
        elif self.settings['model'] == 'exponential': #lethals+beneficials+deleterious
            for i, j in zip(to_fill[0], to_fill[1]):
                randnr = random.random()
                if randnr < self.settings['parameters']['fl']:
                    fitness[i,j] = 0
                elif randnr < self.settings['parameters']['fl']+self.settings['parameters']['fn']:
                    fitness[i,j] = 1
                elif randnr < self.settings['parameters']['fl']+self.settings['parameters']['fn']+self.settings['parameters']['fb']:
                    fitness[i,j] = 1+np.random.exponential(self.settings['parameters']['lb'])
                else:
                    fitness[i,j] = 1-np.random.exponential(self.settings['parameters']['ld'])
                if fitness[i,j] < 0:
                    fitness[i,j] = 0

        elif self.settings['model']  == 'spikes':
            n_spikes = (len(self.settings['parameters']['loc']))

            for i, j in zip(to_fill[0], to_fill[1]):
                randnr = random.random()
                prob = 0
                for spike in range(n_spikes):
                    prob += self.settings['parameters']['freq'][spike]
                    if randnr < prob:
                        fitness[i, j] = self.settings['parameters']['loc'][spike]
                        break
        elif self.settings['model']  == 'lognormal':
            for i, j in zip(to_fill[0], to_fill[1]):
                randnr = random.random()
                if randnr > self.settings['parameters']['fl']:
                    fitness[i, j] = np.random.lognormal(self.settings['parameters']['mu'],
                                                        self.settings['parameters']['sigma'])
        elif self.settings['model']  == 'beta':
            randnr = random.random()
            if randnr > self.settings['parameters']['fl']:
                fitness[i, j] = np.random.lognormal(self.settings['parameters']['a'],
                                                    self.settings['parameters']['b'])
        elif self.settings['model'] == 'from_data':
            for i, j in zip(to_fill[0], to_fill[1]):
                index = np.random.randint(len(self.settings['parameters']['values']))
                if self.settings['parameters']['SD'][index] == 0:
                    new_number = self.settings['parameters']['values'][index]
                else:
                    new_number = np.random.normal(self.settings['parameters']['values'][index],
                                                  self.settings['parameters']['SD'][index])
                if new_number >= 0:
                    fitness[i,j] = new_number
                else:
                    fitness[i,j] = 0

        return fitness

    def get_fitness_effect(self, location, target_base):
        return self.fitness_table[target_base,location]

    def get_nr_offspring(self, sequence_id, return_fitness=False):
        """returns the number of offspring of a sequence according to the fitness
        of that sequence"""
        R0 = self.settings['R0']
        changes = self.current_gen.get_seq(sequence_id)
        fitness = 1
        if changes is not None:
            for pos, base in zip(changes[:, 0], changes[:, 1]):
                fitness *= (self.fitness_table[int(base), int(pos)])
        if return_fitness:
            return np.random.poisson(fitness*R0), fitness
        return np.random.poisson(fitness*R0)

    def mutate_seq(self, pop, seq_id):
        """mutates a sequence (with existing mutations) of length N, according to
        the per base mutation rate"""
        try:
            nr_mut = self.mutations_per_seq.next()
        except StopIteration:
            self.mutations_per_seq = self.new_mutations_per_seq()
            nr_mut = self.mutations_per_seq.next()
        if nr_mut>0:
            success_mut = 0

            while success_mut < nr_mut: #do the mutations one by one
                where = random.randrange(0, self.sequence.len) #draw where the mutation will take place
                base = self.current_gen.get_base(seq_id, where)
                rand_nr = random.random() #draw a random nr for base substitution

                to_check = self.settings['subs_matrix'][int(base), :] #get the cumulative distribution
                                                            #of base substitutions
                new_base = np.where(to_check > rand_nr)[0][0] #find the new base
                if base != new_base:

                    if (base == 1) and (new_base == 0): #G-A mutations
                        if (self.settings['ga_increase']> 1) or (random.random() < self.settings['ga_increase']):
                            pop.add_change(seq_id, where, new_base)
                            success_mut += 1
                    elif (base != 1) or (new_base != 0): #not G-A mutation
                        if (self.settings['ga_increase'] < 1) or random.random() < (1.0/self.settings['ga_increase'] ):
                            pop.add_change(seq_id, where, new_base)
                            success_mut += 1
#    @profile
    def new_generation(self):
        """create a new generation in the simulation"""
        self.effective_pop = 0
        self.gen += 1
        new_generation = Population(self, n_seq=0)
        all_offspring = []
        fitnesses = [0]*self.current_gen.n_seq
        weights = [0]*self.current_gen.n_seq
        #generate offspring list
        for i in xrange(self.current_gen.n_seq):
            #find changes in current sequence
            #find the number of offspring based on the mutations that already took place
            n_offspring, fitness = self.get_nr_offspring(i, return_fitness=True)
            weights[i] = n_offspring
            fitnesses[i] = fitness
         #   temp = [i]*(n_offspring)


            #all_offspring += temp #make a list off the offspring per sequence
        #get average fitness of this generation
        self.average_fitness = np.mean(fitnesses)


        #reduce population to max_pop
        #if len(all_offspring) > self.max_pop:

        if sum(weights) > self.settings['max_pop']:
            #reduce the population randomly to max_pop
            all_offspring = sorted(np.random.choice(range(self.current_gen.n_seq),
                                                    size=int(self.settings['max_pop']),
                                                    p=np.array(weights,dtype=float)/sum(weights)))
        else:
            all_offspring = [i for i,j in enumerate(weights) for k in xrange(j)]
        #actually create the next generation
        ancestor = -1
        for i in all_offspring:
            if i != ancestor:
                ancestor = i
                self.effective_pop += 1
                changes = self.current_gen.get_seq(i)
            seq_id = new_generation.add_sequence(changes)

            self.mutate_seq(new_generation, seq_id)


        if new_generation.n_seq == 0:
            print 'died out'
            n_gen = self.gen
            self = self.copy(self.settings['name'])
            for i in range(n_gen):
                self.new_generation()

        else:
            self.current_gen = new_generation



    def copy(self, name, n_seq = 0):
        '''create a deep copy of the simulation'''
        if n_seq == 0: #original state
            return Simulation(self.settings, sequence = self.sequence,
                              fitness_table=self.fitness_table,name=name)
        else:
            simcopy = Simulation(self.settings, sequence = self.sequence,
                                 fitness_table=self.fitness_table,
                                 name= name, n_seq_init=n_seq)
            sample = self.current_gen.get_sample(n_seq)
            for i,s in enumerate(sample):
                changes = self.current_gen.get_seq(s)
                if changes is not None:
                    for c in changes:
                        simcopy.current_gen.add_change(i,c[0],c[1])
            return simcopy



class Population():
    """class representing a population beloning to a simulation """
    def __init__(self, simulation, changes=None, changed=None, n_seq=1):
        if changes is None:
            self.changed = set([])
            self.changes = {}
            #seqID, location of mutation, mutated target base, ancestor seq ID
        else:
            self.changes = changes
            self.changed = set(changed)
        self.sim = simulation
        self.n_seq = n_seq

    def __len__(self):
        return self.n_seq

    def __str__(self):
        string = '#mutID (from-pos-to)\tsequence\tpatient\n'
        for i in range(self.n_seq):
            if i in self.changed:
                for j in self.changes[i]:
                    pos = j[0]
                    string += '{orig}-{pos}-{to}\t{seq}\t{patient}\n'.format(orig=self.sim.sequence[pos],
                                                                             pos=pos,
                                                                             to=j[1],
                                                                             seq=i,
                                                                             patient=self.sim.settings['name'])
        return string

    def copy(self):
        return Population(self.sim,self.changes,self.changed,self.n_seq)

    def print_sample(self, seq_ids):
        string = '#mutID (from-pos-to)\tsequence\tpatient\n'
        for i in range(self.n_seq):
            if i in self.changed and i in seq_ids:
                for j in self.changes[i]:
                    pos = j[0]
                    string += '{orig}-{pos}-{to}\t{seq}\t{patient}\n'.format(orig=self.sim.sequence[pos],
                                                                             pos=pos,
                                                                             to=j[1],
                                                                             seq=i,
                                                                             patient=self.sim.settings['name'])
        print string

    def get_sample(self, sample_size):
        try:
            return np.random.choice(self.n_seq,size=sample_size,replace=False)
        except ValueError:
            return range(self.n_seq)




    def add_sequence(self, changes=None):
        self.n_seq += 1
        if changes is not None:
            self.changed.add(self.n_seq-1)
            self.changes[self.n_seq-1] = changes

        return self.n_seq-1

    def add_change(self, seq_id, pos, target):
        if seq_id in self.changed:
            #add to existing changes list
            if pos in self.changes[seq_id][:, 0]:
                self.changes[seq_id][self.changes[seq_id][:, 0] == pos, 1] = target
            else:
                self.changes[seq_id] = np.vstack((self.changes[seq_id], [pos, target]))
        else:
            #add a new changed sequence
            self.changed.add(seq_id)
            self.changes[seq_id] = np.array([[pos, target]])

    def get_base(self, seq_id, pos):
        if seq_id in self.changed:
            if pos in self.changes[seq_id][:, 0]:
                return self.changes[seq_id][self.changes[seq_id][:, 0] == pos, 1]
        return self.sim.sequence[pos]

    def stats(self):
        stats = {}
        stats['n_seq'] = self.n_seq
        stats['unmutated'] = self.n_seq-len(self.changes)
        if len(self.changed)>0:
            all_mutations = np.vstack(self.changes.values())
        else:
            all_mutations = []
        stats['total_mutations'] = len(all_mutations)
        all_mutations = [tuple(row) for row in all_mutations]
        stats['unique_mutations'] = len(set(all_mutations))

        mut_counts = np.array(Counter(all_mutations).values())
        if len(mut_counts) > 0:
            stats['majority_mutations'] = sum(mut_counts > (stats['n_seq']/2.0))
            stats['max_fraction'] = max(mut_counts/float(stats['n_seq']))
        else:
            stats['majority_mutations'] = 0
            stats['max_fraction'] = 0
        GA = 0
        for i in all_mutations:
            if self.sim.sequence[i[0]]==1 and i[1] == 0:
                GA+=1.0
        try:
            stats['GA_rate'] = GA/len(all_mutations)
        except ZeroDivisionError:
            stats['GA_rate'] = None
        return stats

    def to_fasta(self, seq_ids=[], n_seq=1, description=''):
        string = ''
        if len(seq_ids) == 0:
            seq_ids = random.sample(range(self.n_seq), n_seq)
        bar = progressbar.ProgressBar()
        for i in bar(range(len(seq_ids))):
            seqID = seq_ids[i]
            string += '>'+str(seqID)+''+str(description)+'\n'
            changed_here = self.get_seq(seqID)
            seq = deepcopy(self.sim.sequence)
            if changed_here is not None:
                for i in changed_here:
                    seq.sequence[int(i[0])] = int(i[1])
            string += str(seq)+'\n'
        return string

    def consensus_sequence(self):
        seq = deepcopy(self.sim.sequence)
        all_mutations = np.vstack(self.changes.values())
        all_mutations = [tuple(row) for row in all_mutations]

        mutations = Counter(all_mutations)
        for mut in mutations.keys():
            if mutations[mut] >= self.n_seq/2.0:
                seq.sequence[int(mut[0])] = int(mut[1])
        return seq

    def get_seq(self, sequence_id):
        if sequence_id in self.changed:
            return self.changes[sequence_id]
        else:
            return None


if __name__ == '__main__':
    import sys

    settings = sys.argv[1]
    n_gen = int(sys.argv[2])

    sim = Simulation(settings)

    for i in range(n_gen):
        sim.new_generation()

    print sim
