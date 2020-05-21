"""
Created on Thu Dec 15 13:42:33 2016

@author: eva
"""
from __future__ import print_function
import pkg_resources
import sys
from collections import Counter
import random
from copy import copy, deepcopy
import time
import os
import inspect

import numpy as np
np.set_printoptions(threshold=sys.maxsize, suppress=True)
import yaml
import scipy.stats as scats

class Seq(object):
    '''
    Sequence class.

    **Attributes**:
    * `translation` (str): the translation order from bases (ATGC) to numbers
    * `len` (int): the length of the sequence
    * `sequence` (list): a list representation of the sequence

    **Keyword Arguments**:
    * `seq_len` (int): length of the sequence to generate. Will be ignored if `seq` is provided. Default 100.
    * `base_dist` (list): list of length 4, base distribution in order A-G-T-C.
    Distribution is cumulative: equal base distribution would be represented
    as [0.25,0.5,0.75,1]. When not provided, will default to equal base distribution
    * `seq` (str): the sequence to use. When not provided, a random sequence
    will be generated according to `seq_len` and `base_dist`
    '''
    def __init__(self, seq_len=100, base_dist=None, seq='',**kwargs):
        '''
        Create a sequence object.

        *Keyword Arguments*:
        * `seq_len` (int): length of the sequence to generate. Will be ignored if `seq` is provided. Default 100.
        * `base_dist` (list): list of length 4, base distribution in order A-G-T-C.
        Distribution is cumulative: equal base distribution would be represented
        as [0.25,0.5,0.75,1]. When not provided, will default to equal base distribution
        * `seq` (str): the sequence to use. When not provided, a random sequence
        will be generated according to `seq_len` and `base_dist`  '''
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
        return self.translation[self.sequence[key]]

    def translate_sequence(self, seq):
        '''
        translate a sequence from bases to numbers

        **Arguments**:
        * `seq`: string (base) representation of the sequence

        **Returns**:
        list (number) representation of the sequence
        '''
        sequence = []
        for i in seq:
            sequence.append(self.translation.index(i))
        return sequence

    def generate_seq(self, seq_len, base_dist=None):
        '''
        Sets the sequence to a random sequence of seq_len bases according to the
        base distribution

        **Arguments**:
        * `seq_len` (int): the length of the sequence to generate
        * `base_dist` (list): list of length 4, base distribution in order A-G-T-C.
        Distribution is cumulative: equal base distribution would be
        represented as [0.25,0.5,0.75,1] (which is the default)

        **Returns**:
            Nothing.
        '''
        if seq_len < 0:
            raise TypeError('seq_len must be positive integer')
        seq = []
        if base_dist is None:
            base_dist = np.array([0.25, 0.5, 0.75, 1])
        else:
            base_dist = np.array(base_dist)
        for i in range(seq_len):
            seq.append(4-sum(base_dist > random.random()))
        self.sequence = seq

class Simulation(object):
    '''
    Sequence simulation object

    **Attributes**:
    * `settings` (Dict): all simulation settings
    * `mut_rate` (float): mutation rate for the simulation (adapted for potential
    increase in G-A rate)
    * `sequence` (Seq): ancestor sequence the simulation starts with
    * `fitness_table` (4xseq_len list): table containing fitness values of all
    possble bases in the sequence
    * `gen` (int): the current generation
    * `average_fitness` (float): average fitness of the population
    * `effective_pop` (int): effective population size of current generation
    (number of unique sequences)
    * `n_seq` (int) total number of sequenes in the current generation
    * `current_gen` (Population): the current generation
    * `mutations_per_seq` (list): number of mutations in the next 1000 mutation
    events. implemented for efficiency.

    **Arguments**:
    * `simulation_settings` (str): either one of the default simulation
    settings (HIV (default), phix174), or a path to a settings file
    in yaml-format

    **note**: loaded simulation settings will be overwritten by parameters
        provided as keyword arguments

    **Keyword Arguments**:
    * `model` (str): model to use for the MFED (otions: 'neutral'(default),
    'lognormal','exponential','spikes','beta', 'from_data','steps')
    * `parameters` (dict): parameters for the MFED
    * `mut_rate` (float): mutation rate per site per generation (default = 2.16e-5)
    * `subs_matrix` (4x4 list): cumulative substitution matrix
    rows = from, columns = to
    * `seq_len` (int): length of sequence to simulate (default 2600).
    used only if no sequence is provided
    * `basedist` (4x1 list): list of cumulutive distribution of bases used
    for generating the sequence in order [A,G,T,C]. default: [0.25,0.5,0.75,1]
    * `R0` (float): initial average amount of offspring per sequence
    * `ga_increase` (float): Increase in G-A mutation rate (defaul: 1 (no increase))
    * `max_pop` (int): maximum population size
    * `name` (str): name of the simulation (used in output)
    * `sequence` (Seq): the sequence to use for the simulation
    * `fitness_table` (4xseq_len): fitness table to use for the simulation
    * `n_seq_init` (int): initial number of sequences (default 1)

    '''

    def __init__(self, simulation_settings='HIV', **kwargs):
        '''
        create a Simulation object.

        Arguments:
            simulation_settings (str): either one of the default simulation
                settings (HIV (default), phix174), or a path to a settings file
                in yaml-format

        **note**: loaded simulation settings will be overwritten by parameters
            provided as keyword arguments

        Keyword Arguments:
            model (str): model to use for the MFED (otions: 'neutral'(default),
                'lognormal','exponential','spikes','beta', 'from_data','steps',
                'lognormal_ben_only')
            parameters (dict): parameters for the MFED
            mut_rate (float): mutation rate per site per generation (default = 2.16e-5)
            subs_matrix (4x4 list): cumulative substitution matrix
                rows = from, columns = to
            seq_len (int): length of sequence to simulate (default 2600).
                used only if no sequence is provided
            basedist (4x1 list): list of cumulutive distribution of bases used
                for generating the sequence in order [A,G,T,C]. default: [0.25,0.5,0.75,1]
            R0 (float): initial average amount of offspring per sequence
            ga_increase (float): Increase in G-A mutation rate (defaul: 1 (no increase))
            max_pop (int): maximum population size
            name (str): name of the simulation (used in output)
            sequence (Seq): the sequence to use for the simulation
            fitness_table (4xseq_len): fitness table to use for the simulation
            n_seq_init (int): initial number of sequences (default 1)
        '''

        if type(simulation_settings) is list:
            simulation_settings = simulation_settings[0]
        try:
            with open(simulation_settings) as f:
                self.settings = yaml.safe_load(f)
        except IOError:
            path = pkg_resources.resource_filename('SeqSimEvo','simulation_settings/')
            with open(path+simulation_settings) as f:
                self.settings = yaml.safe_load(f)
        except TypeError:
            self.settings = simulation_settings

        for key, value in kwargs.items():
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
            self.fitness_table = self.getfitness_table()
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
        '''select the number of mutations that will happen for the next 1000
        mutated sequences. implemented for increased efficiency'''
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


    def getfitness_table(self):
        '''creates a table with random fitness values according to model and
        parameters par in settings'''
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
        elif self.settings['model']  in ['lognormal','lognormal_ben_only']:
            for i, j in zip(to_fill[0], to_fill[1]):
                randnr = random.random()
                if randnr > self.settings['parameters']['fl']:
                    fitness[i, j] = np.random.lognormal(self.settings['parameters']['mu'],
                                                        self.settings['parameters']['sigma'])
            if self.settings['model'] == 'lognormal_ben_only':
                fitness[fitness < 1] = 1
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
        '''
        get the fitness effect of a mutation at location to target_base

        **Arguments**:
        * `location` (int): the position in the genome (0-based counting)
        * `target_base` (int): the base of which to get the fitness effect

        **Returns**:
            the fitness effect at this position (float)
        '''
        return self.fitness_table[target_base,location]

    def get_nr_offspring(self, sequence_id, return_fitness=False):
        """
        get the number of offspring of a sequence according to the fitness
        of that sequence

        **Arguments**:
        * `sequence_id` (int): the id of the sequence to get offspring of
        * `return_fitness` (Bool): return fitness as well

        **Returns**:
        the calculated number of offspring (int) (if return_fitness is False)
        tuple (int,float): calculated number of offsprint and the fitness of
        the sequence (if return_fitness is True)
        """
        R0 = self.settings['R0']
        changes = self.current_gen.get_seq(sequence_id)
        fitness = 1

        #calculate the fitness of the sequence
        if changes is not None:
            for pos, base in zip(changes[:, 0], changes[:, 1]):
                fitness *= (self.fitness_table[int(base), int(pos)])

        if self.settings['offspring_distribution'] == 'poisson':
            offspring = np.random.poisson(fitness*R0)
        elif self.settings['offspring_distribution'] == 'normal':
            offspring = int(round(np.random.normal(loc=fitness*R0,scale=self.settings['offspring_sigma'])))
        else:
            raise ValueError('offspring distribution {} not understood'.format(self.settings['offspring_distribution']))
        if return_fitness:
            return offspring, fitness
        return offspring

    def mutate_seq(self, pop, seq_id_new,seq_id_old):
        """
        mutates a sequence (with existing mutations)

        **Arguments**:
        * `pop` (Population): the population which holds the new generation
        * `seq_id_new` (int): the sequence ID of the new sequence
        * `seq_id_old` (int): the sequence ID of the old sequence
        (in the current generation)

        **Returns**:
            True/False (did the sequence mutate or not?)
        """
        #get the number of mutations that will take place
        try:
            nr_mut = next(self.mutations_per_seq)
        except StopIteration:
            self.mutations_per_seq = self.new_mutations_per_seq()
            nr_mut = next(self.mutations_per_seq)


        if nr_mut>0:
            success_mut = 0
            while success_mut < nr_mut: #do the mutations one by one
                where = random.randrange(0, self.sequence.len) #draw where the mutation will take place
                base = self.current_gen.get_base(seq_id_old, where)
                rand_nr = random.random() #draw a random nr for base substitution

                to_check = self.settings['subs_matrix'][int(base), :]
                    #get the cumulative distribution of base substitutions
                new_base = np.where(to_check > rand_nr)[0][0] #find the new base

                #only accept the mutation if there was an actual change and make
                #sure mutations are accepted in line with the G-A increase
                if base != new_base:
                    if (base == 1) and (new_base == 0): #G-A mutations
                        if (self.settings['ga_increase']>= 1) or (random.random() < self.settings['ga_increase']):
                            pop.add_change(seq_id_new, where, new_base)
                            success_mut += 1
                    elif (base != 1) or (new_base != 0): #not G-A mutation
                        if (self.settings['ga_increase'] <= 1) or random.random() < (1.0/self.settings['ga_increase'] ):
                            pop.add_change(seq_id_new, where, new_base)
                            success_mut += 1
            return True
        else:
            return False

    def new_generation(self,new_gen=None,dieout=False):
        """
        Create a new generation in the simulation

        **Arguments**:
        * `new_gen` (Population): the population to store the new generation in
        (a new Population instance will be created if left out)
        * `dieout` (Bool): allow the populations to die out (True) or redo
        new_generation if the population died out (False, default)

        **Returns**:
            Nothing. Updates current_gen, effective_pop, gen, average_fitness,
            and n_seq
        """
        self.effective_pop = 0
        self.gen += 1
        if new_gen is None:
            new_gen = Population(self, n_seq=0)
        all_offspring = []

        fitnesses = [0]*self.current_gen.n_seq #will hold the fitness values for all sequences
        weights = [0]*self.current_gen.n_seq #will hold the weights (number of ofspring) for all sequences

        #generate offspring list
        try:
            xrange()
        except NameError:
            xrange = range

        for i in xrange(self.current_gen.n_seq):
            #find changes in current sequence
            #find the number of offspring based on the mutations that already took place
            n_offspring, fitness = self.get_nr_offspring(i, return_fitness=True)
            weights[i] = n_offspring
            fitnesses[i] = fitness

        #get average fitness of this generation
        self.average_fitness = np.mean(fitnesses)

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
            seq_id = new_gen.add_sequence(changes=changes)

            self.mutate_seq(new_gen, seq_id, i)


        if new_gen.n_seq == 0 and not dieout :
            print('died out')
            n_gen = self.gen
            self = self.copy(self.settings['name'])
            for i in range(n_gen):
                self.new_generation()

        else:
            self.current_gen = new_gen
        self.n_seq = self.current_gen.n_seq



    def copy(self, name, n_seq = -1,**kwargs):
        '''
        create a deep copy of the simulation, number of generations will be set
        back to 0

        **Arguments**:
        * `name`: the name to use for the copy
        * `n_seq`: the number of sequences to keep in the copy. -1: original state

        **Returns**:
            a copy of the simulation in its original state
        '''
        if n_seq == -1: #original state
            return Simulation(deepcopy(self.settings), sequence = self.sequence,
                              fitness_table=self.fitness_table,name=name,**kwargs)
        else:
            simcopy = Simulation(deepcopy(self.settings), sequence = self.sequence,
                                 fitness_table=self.fitness_table,
                                 name= name, n_seq_init=n_seq,**kwargs)
            sample = self.current_gen.get_sample(n_seq)
            for i,s in enumerate(sample):
                changes = self.current_gen.get_seq(s)
                if changes is not None:
                    for c in changes:
                        simcopy.current_gen.add_change(i,c[0],c[1])
            return simcopy



class Population():
    """
    class representing a population beloning to a simulation

    **Attributes**:
    * `changed` (set): the sequence IDs that countain a mutation
    * `changes` (dict): contains a numpy array per sequence with the position
    and new base per mutation
    * `sim` (Simulation): the simulation to which this Population belongs
    * `n_seq`: the number of sequences in the population

    **Arguments**:
    * `simulation` (Simulation): the simulation to which this Population belongs
    * `changes` (dict): for pre-exsisting changes. Numpy array per sequence with
    the position and new base per mutation
    * `changed` (set): for pre-existing changes: the sequence IDs that countain a mutation
    * `n_seq` (int): number of sequences in the population (default 1)

    """
    def __init__(self, simulation, changes=None, changed=None, n_seq=None,**kwargs):
        """
        Initialize the Population

        Arguments:
            simulation (Simulation): the simulation to which this Population belongs
            changes (dict): for pre-exsisting changes. Numpy array per sequence with
                the position and new base per mutation
            changed (set): for pre-existing changes: the sequence IDs that countain a mutation
            n_seq (int): number of sequences in the population (default 1)
        """
        if changes is None:
            self.changed = set([])
            self.changes = {}
        else:
            self.changes = changes
            for i in self.changes:
                self.changes[i] = np.array(self.changes[i])
            self.changed = set(changed)
        self.sim = simulation
        if n_seq is None:
            self.n_seq = self.sim.n_seq
        else:
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
                                                                             to=self.sim.sequence.translation[j[1]],
                                                                             seq=i,
                                                                             patient=self.sim.settings['name'])
        return string

    def copy(self):
        ''' create a deep copy of the population'''
        return Population(self.sim.copy(str(self.sim.settings['name'])+'_copy'),deepcopy(self.changes),deepcopy(self.changed),self.n_seq)

    def print_sample(self, seq_ids):
        ''' print a summary of the mutations that have occured

        Prints a summary of the mutations that have occured in all seq_ids in
        the format: #mutID (from-pos-to)\tsequence\tpatient\n

        **Arguments**:
        * `seq_ids` (list): the ids of the sequences to print a summary of

        **Returns**:
        Nothing. Output is printed to stdout.
        '''
        if any(np.array(seq_ids)>self.n_seq):
            raise IndexError('seqID out of range')
        string = '#mutID (from-pos-to)\tsequence\tpatient\n'
        for i in range(self.n_seq):
            if i in self.changed and i in seq_ids:
                for j in self.changes[i]:
                    pos = j[0]
                    string += '{orig}-{pos}-{to}\t{seq}\t{patient}\n'.format(orig=self.sim.sequence[pos],
                                                                             pos=pos,
                                                                             to=self.sim.sequence.translation[j[1]],
                                                                             seq=i,
                                                                             patient=self.sim.settings['name'])
        print(string)

    def sample_to_string(self, seq_ids):
        ''' create a summary of the mutations that have occured

        **Arguments**:
        * `seq_ids` (list): the ids of the sequences to print a summary of

        **Returns**:
        * `string`: summary of the mutations that have occured in the seq_ids,
        in the format "#mutID (from-pos-to)\tsequence\tpatient"

        '''
        string = '#mutID (from-pos-to)\tsequence\tpatient\n'
        for i in range(self.n_seq):
            if i in self.changed and i in seq_ids:
                for j in self.changes[i]:
                    pos = j[0]
                    string += '{orig}-{pos}-{to}\t{seq}\t{patient}\n'.format(orig=self.sim.sequence[pos],
                                                                             pos=pos,
                                                                             to=self.sim.sequence.translation[j[1]],
                                                                             seq=i,
                                                                             patient=self.sim.settings['name'])
        return string

    def get_sample(self, sample_size):
        ''' get a random sample from the population
        If the sample size is larger than the population, the whole population is
        returned

        **Arguments**:
        * `sample_size` (int): the size of the sample

        **Returns**:
            a list of sequence IDs randomly sampled from the population
        '''
        try:
            return np.random.choice(self.n_seq,size=int(sample_size),replace=False)
        except ValueError:
            return range(self.n_seq)


    def delete_sequence(self, ID):
        ''' delete sequence from the population.
        Sequence IDs will be reassigned to fit within the new number of sequences.

        **Arguments**:
        * `ID` (int): the sequence ID of the sequence to remove

        **Returns**:
        nothing. Sequence is removed in-place.
        '''
        self.n_seq-=1
        if self.get_seq(ID) is not None:
            self.changed.remove(ID)
            del self.changes[ID]

        if self.get_seq(self.n_seq) is not None:
            self.changed.remove(self.n_seq)
            self.changed.add(ID)

            self.changes[ID] = self.changes[self.n_seq]
            del self.changes[self.n_seq]


    def add_sequence(self, changes=None):
        '''add a sequence to the population

        add a sequence, optionally with certain changes (as a list of position, new), to the population

        **Arguments**:
        * `changes`: the changes present in this sequence

        **Returns**:
        the sequence ID of the newly added sequence
        '''
        self.n_seq += 1
        if changes is not None:
            for i in changes:
                self.add_change(self.n_seq-1,i[0],i[1])

        return self.n_seq-1

    def add_change(self, seq_id, pos, target):
        ''' add a change to an existing sequence

        **Arguments**:
        * `seq_id` (int): the sequence ID to add the change to
        * `pos` (int): the position of the change
        * `target` (int): the new base at the changed position

        **Returns**:
            Nothing. Population is changed in-place.
        '''
        if pos > len(self.sim.sequence):
            raise IndexError('Pos {} outside sequence length'.format(pos))

        if seq_id > self.n_seq:
            raise IndexError('SeqID {} outside pop size {} {}'.format(seq_id, self.n_seq,self))


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
        ''' get the current base at position pos in sequence with id seq_id '''
        if seq_id in self.changed:
            if pos in self.changes[seq_id][:, 0]:
                return self.changes[seq_id][self.changes[seq_id][:, 0] == pos, 1]
        return self.sim.sequence.sequence[pos]

    def stats(self):
        ''' return a dict of stats about the population

        **Keys in the returned dict:**
        * `n_seq`: the total number of sequences in the current generation
        * `unmutated`: the number of unmutated sequences
        * `total_mutations`: the number of mutations in total
        * `unique_mutations`: the length of the set of all mutations
        * `majority_mutations`: the number of mutations that reached majority
        * `max_fraction`: the highest fraction reached by a mutation
        * `GA_rate`: the fraction of mutations that are G-to-A
        '''
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

        if len(all_mutations) > 0:
            mut_counts = np.array(list(Counter(all_mutations).values()))
        else:
            mut_counts = []
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


    def to_fasta(self, seq_ids=[], n_seq=None, description='',progress=False):
        '''
        convert (part of) the population to fasta-format.

        Without any arguments, all sequences in the population will be returned.

        **Optional arguments**:
        * `seq_ids` (list): list of sequence IDs to convert to fasta-format
        * `n_seq`: number of sequences to convert to fasta-format (random draw
        from population). This number will be ignored if seq_ids is given
        * `description` (str): description of the sequences, will be added after
        the sequenceID in the header of each sequence
        * `progress` (Bool): display a progress bar?

        **Returns**:
            str: the selected sequences in fasta-format
        '''
        string = ''
        if len(seq_ids) == 0:
            if n_seq is None or n_seq > self.n_seq:
                n_seq = self.n_seq

            seq_ids = random.sample(range(self.n_seq), n_seq)
        for i in range(len(seq_ids)):
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
        ''' return the consensus sequence of the population'''
        seq = deepcopy(self.sim.sequence)
        all_mutations = np.vstack(self.changes.values())
        all_mutations = [tuple(row) for row in all_mutations]

        mutations = Counter(all_mutations)
        for mut in mutations.keys():
            if mutations[mut] >= self.n_seq/2.0:
                seq.sequence[int(mut[0])] = int(mut[1])
        return seq

    def get_seq(self, sequence_id):
        ''' get the changes in the sequence with id sequence_id

        **Raises**:
            `IndexError`: when sequence_id is out of bounds
        '''
        if sequence_id > self.n_seq:
            raise IndexError('sequence_id is out of bounds')
        elif sequence_id in self.changed:
            return self.changes[sequence_id]
        else:
            return None

    def Hamming_distance(self,sample,action='mean'):
        '''
        calculate the inter-sequence hamming distances in a sample.
        if action is 'mean', return the mean hamming distance,
        if action is 'Poisson_fit', return the poisson fit for the time since
            infection from the distribution of hamming distances as presented
            in Lee et al, 2010
        '''
        simulation_settings = self.sim.settings
        HDs = []
        for i in sample:
            if i in self.changed:
                changed1 = [str(k) for k in self.changes[i]]
            else:
                changed1 = []
            for j in sample:
                if i!=j:
                    if j in self.changed:
                        changed2 =  [str(k) for k in self.changes[j]]
                    else: changed2 = []
                    HDs.append(len(set(list(changed1)) ^ set(list(changed2))))
        if action == 'mean':
            return np.mean(HDs)
        elif action == 'Poisson_fit':
            poiss = np.mean(HDs)/(2*simulation_settings['mut_rate']*simulation_settings['seq_len'])
            exp = scats.poisson.pmf(range(max(HDs)+1),np.mean(HDs))*len(HDs)
            obs = np.histogram(HDs, bins=range(0,max(HDs)+2))[0]
            pval =scats.chisquare(obs,exp,ddof=len(exp)-1-len(sample)).pvalue
            if np.isnan(pval) or pval>0.05:
                return poiss
            else:
                return np.nan
