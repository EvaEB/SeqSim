from __future__ import print_function
import SeqSimEvo as seq_sim
import numpy as np
from collections import Counter
import matplotlib.pyplot as plt
import random
import scipy.stats as scats
import os

class Simulation(seq_sim.Simulation):
    def __init__(self,simulation_settings='HIV',tag_dist=None,tag_len=8,n_seq_init=1,**kwargs):
        seq_sim.Simulation.__init__(self,simulation_settings=simulation_settings,
                                    n_seq_init=n_seq_init,**kwargs)

        self.sequence = Seq(tag_dist,tag_len,seq_len=self.settings['seq_len'],
                            base_dist=self.settings['basedist'],**kwargs)
        self.current_gen = Population(self,n_seq=int(self.settings['n_seq_init']),**kwargs)
        self.mutations_per_tag = self.new_mutations_per_tag()


    def new_mutations_per_tag(self):
        '''select the number of mutations that will happen in tags for the next 1000
        mutated sequences. implemented for increased efficiency'''
        return iter(np.random.binomial(self.sequence.tag_len,self.settings['mut_rate'],1000))


    def mutate_seq(self,pop,seq_id_new,seq_id_old):
        #do the original mutation as before
        seq_sim.Simulation.mutate_seq(self,pop, seq_id_new,seq_id_old)
        #print self.current_gen.get_seq(seq_id_old), pop.get_seq(seq_id_new)
        tag = self.current_gen.tags[seq_id_old]
        #add new functionality: mutate tag
        try:
            nr_mut = next(self.mutations_per_tag)
        except StopIteration:
            self.mutations_per_tag = self.new_mutations_per_tag()
            nr_mut = next(self.mutations_per_tag)
        if nr_mut>0:
            tag_seq = self.sequence.get_tag(tag)
            success_mut = 0
            while success_mut < nr_mut: #do the mutations one by one
                where = random.randrange(0, self.sequence.tag_len) #draw where the mutation will take place

                base = tag_seq[where]
                rand_nr = random.random() #draw a random nr for base substitution

                to_check = self.settings['subs_matrix'][int(base), :] #get the cumulative distribution
                                                                      #of base substitutions
                new_base = np.where(to_check > rand_nr)[0][0] #find the new base
                if base != new_base:

                    if (base == 1) and (new_base == 0): #G-A mutations
                        if (self.settings['ga_increase']> 1) or (random.random() < self.settings['ga_increase']):
                            tag_seq[where] = new_base
                            new_tag = self.sequence.to_tag(tag_seq)
                            pop.tags[seq_id_new] = new_tag
                            success_mut += 1
                    elif (base != 1) or (new_base != 0): #not G-A mutation
                        if (self.settings['ga_increase'] < 1) or random.random() < (1.0/self.settings['ga_increase'] ):
                            tag_seq[where] = new_base
                            new_tag = self.sequence.to_tag(tag_seq)
                            success_mut += 1
        else:
            new_tag = tag
        pop.tags[seq_id_new] = new_tag



class Population(seq_sim.Population):
    def __init__(self,simulation,tags=None,**kwargs):
        seq_sim.Population.__init__(self,simulation,**kwargs)
        if tags is None:
            pos_tags = np.random.choice(range(4**self.sim.sequence.tag_len),len(self.sim.sequence.tag_dist),replace=False)
            self.tags = np.random.choice(pos_tags,self.n_seq,p=self.sim.sequence.tag_dist)
        else:
            self.tags = tags

    def get_tag(self,seqID):
        return self.sim.sequence.get_tag(self.tags[seqID])

    def add_sequence(self,tag=-1,changes=None):
        #print changes
        new_seq_id = seq_sim.Population.add_sequence(self,changes=changes)
        self.tags.append(tag)
        return new_seq_id

    def print_sample(self,seq_ids):
        ''' print a summary of the mutation that have occured in all seq_ids in
        the format: #mutID (from-pos-to)\tsequence\tpatient\n'''
        string = '#mutID (from-pos-to)\tsequence\tpatient\n'
        for i in range(self.n_seq):
            if i in self.changed and i in seq_ids:
                for j in self.changes[i]:
                    pos = j[0]
                    string += '{orig}-{pos}-{to}\t{seq}\t{patient}\t{tag}\n'.format(orig=self.sim.sequence[pos],
                                                                                    pos=pos,
                                                                                    to=j[1],
                                                                                    seq=i,
                                                                                    patient=self.sim.settings['name'],
                                                                                    tag = self.tags[i])
        print(string)



class Seq(seq_sim.Seq):
    def __init__(self,tag_dist=None,tag_len=8,**kwargs):
        seq_sim.Seq.__init__(self,**kwargs)
        if tag_dist is None:
            self.tag_dist = np.ones(4**tag_len)/(4.0**tag_len)
        else:
            self.tag_dist = tag_dist
        self.tag_len = tag_len

    def get_tag(self,number,to_return='list'):
        '''get the translation of the tag into a list of numbers
        (``return = 'list'``) or as a string (``return = 'seq'``)'''
        if to_return == 'seq':
            bases = 'ACGT'
            tag = ''
        else:
            tag = []
        for i in range(self.tag_len):
            if to_return == 'seq':
                tag =bases[number%4]+tag
            else:
                tag = [number%4] + tag
            number = number/4
        return tag

    def to_tag(self,lst):
        '''get the tag number from a list of bases'''
        number = 0
        power = 0
        for i in lst[::-1]:
            number += (4**power)*i
            power+=1
        return number

def output_tags(sim,transfer):
    print('#transfer', transfer)
    counts = Counter(sim.current_gen.tags)
    for i in counts:
        print(i, counts[i])

def output_fitness_per_tag(sim,transfer):
    #print '#',transfer
    tags = {}
    # changes = []
    for i in range(sim.current_gen.n_seq):
        fitness =  sim.get_nr_offspring(i,return_fitness=True)[1]
        tag = sim.current_gen.tags[i]
        try:
            tags[tag].append(fitness)
        except KeyError:
            tags[tag] = [fitness]
        # try:
        #     changes.append(len(sim.current_gen.get_seq(i)))
        # except TypeError:
        #     changes.append(0)

    for i in tags:
        print(transfer, i, np.mean(tags[i]), np.max(tags[i]),np.min(tags[i]),len(tags[i]))

def output_seqs(sim,transfer):
    if transfer in [99,199]:
        print('#',transfer)
        sim.current_gen.print_sample(sim.current_gen.get_sample(sim.current_gen.n_seq))


def run_sim(tag_dist, pop_size,transfer_prop,n_transfer,div,scenario,output,progress=False,param=None):
    if param is None:
        sim = Simulation(simulation_settings='phix174',tag_dist=tag_dist,
                         n_seq_init=pop_size,model=scenario,max_pop=pop_size)
    else:
        sim = Simulation(simulation_settings='phix174',tag_dist=tag_dist,
                         n_seq_init=pop_size,model=scenario,max_pop=pop_size,
                         parameters=param)
    sim.settings['R0'] = 110

    #diversify sequences
    for i in range(sim.current_gen.n_seq):
        #how many changes
        n_change = np.random.poisson(div*sim.sequence.len)
        for change in range(n_change):
            where = int(np.random.uniform(0,sim.sequence.len))
            while True:
                new = int(np.random.uniform(0,4))
                if new != sim.sequence[where]:
                    break
            sim.current_gen.add_change(i, where, new)

    if progress:
        pbar = tqdm.tqdm(total=n_transfer)
    for i in range(n_transfer):
        output(sim,i)
        sim.new_generation(new_gen=Population(sim,tags=[],n_seq=0))
        sim.settings['max_pop'] = pop_size

        sim.new_generation(new_gen=Population(sim,tags=[],n_seq=0))
        sim.settings['max_pop'] =transfer_prop*sim.current_gen.n_seq
        if progress:
            pbar.update(1)
    if progress:
        pbar.close()

    output(sim, i+1)
    return sim

def output_all(sim, transfer):
    count = 0
    changes_old = ''
    tag_old= ''
    for seqID in range(len(sim.current_gen)):
        changes = sim.current_gen.get_seq(seqID)
        if changes is not None:
            changes_string = '/'.join(['{}-{}'.format(i[0],i[1]) for i in changes])
        else:
            changes_string = ''
        if (sim.current_gen.tags[seqID] == tag_old) and (changes_string == changes_old):
            count+=1
        else:
            tag_old = sim.current_gen.tags[seqID]
            changes_old = changes_string
            data = {'count': count,
                    'tag': sim.current_gen.tags[seqID],
                    'fitness': sim.get_nr_offspring(seqID,return_fitness=True)[1],
                    'transfer': transfer,
                    'experimentID':experimentID,
                    'changes': changes_string
                   }
            print('{count},{tag},{fitness},{transfer},{experimentID},{changes}'.format(**data))
            count=1
    data = {'count': count,
            'tag': sim.current_gen.tags[seqID],
            'fitness': sim.get_nr_offspring(seqID,return_fitness=True)[1],
            'transfer': transfer,
            'experimentID':experimentID,
            'changes': changes_string
           }
    print('{count},{tag},{fitness},{transfer},{experimentID},{changes}'.format(**data))


if __name__ == '__main__':
    import sys
    import tqdm

    if sys.argv[1] in ['help','-h','h']:
        print('usage: python sequence_tags.py MFED nBarcodes distBarcodes popSize diverseStart transferProp')
        exit()
    elif sys.argv[1] == 'neutral':
        scenario = 'neutral'
        param=[]
    elif sys.argv[1] == 'selection':
        scenario = 'exponential'
    elif sys.argv[1] == 'more_lethal':
        scenario = 'exponential'
        param = {'fl': 0.4,'fb': 0.2,'lb': 0.03,'fd': 0.4, 'ld': 0.21}
    elif sys.argv[1] == 'more_ben':
        scenario = 'exponential'
        param = {'fl': 0.2,'fb': 0.4,'lb': 0.03,'fd': 0.4, 'ld': 0.21}
    elif sys.argv[1] == 'more_neutral':
        scenario = 'exponential'
        param = {'fl': 0.1,'fb': 0.3,'fn':0.3,'lb': 0.03,'fd': 0.3, 'ld': 0.21}
    elif sys.argv[1] == 'no_beneficial':
            scenario = 'exponential'
            param = {'fl': 0.1,'fb': 0.0,'fn':0.6,'lb': 0,'fd': 0.3, 'ld': 0.21}
    else:
        print('unknown MFED')
        exit()

    n_tags = int(sys.argv[2])

    if sys.argv[3] == 'exponential':
        tag_dist = scats.expon.pdf(range(n_tags),scale=n_tags/3) # -->exponential
    elif sys.argv[3] == 'uniform':
        tag_dist = np.ones(n_tags)
    else:
        raise('unknown tag distribution')
    tag_dist = tag_dist/sum(tag_dist)

    pop_size = int(sys.argv[4])

    div = float(sys.argv[5]) #before setting this != 0: make sure diversified sequences don't have too large fitness effects!!!!!!!
    if div != 0:
        print('diversity diffferent from 0 not implemented yet!!')
        exit()

    transfer_prop = float(sys.argv[6])

    n_transfer = 10

    experimentID = '_'.join(sys.argv[1:])
    #run_sim(tag_dist, pop_size,transfer_prop,n_transfer,div,scenario,output_tags)
    #sim=run_sim(tag_dist, pop_size,transfer_prop,n_transfer,div,scenario,output_fitness_per_tag)
    #run_sim(tag_dist, pop_size,transfer_prop,n_transfer,div,scenario,output_seqs,progress=False)
    sim = run_sim(tag_dist, pop_size, transfer_prop, n_transfer, div, scenario, output_all,param=param)
