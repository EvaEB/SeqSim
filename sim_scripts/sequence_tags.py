import seq_sim
import numpy as np
from collections import Counter
import matplotlib.pyplot as plt
import random
import scipy.stats as scats

class Simulation(seq_sim.Simulation):
    def __init__(self,tag_dist=None,tag_len=8,n_seq_init=1,**kwargs):
        self.sequence = Seq(tag_dist,tag_len,**kwargs)
        seq_sim.Simulation.__init__(self,sequence=self.sequence,n_seq_init=n_seq_init,**kwargs)
        self.current_gen = Population(self,n_seq=int(self.settings['n_seq_init']),**kwargs)
        self.mutations_per_tag = self.new_mutations_per_tag()


    def new_mutations_per_tag(self):
        '''select the number of mutations that will happen in tags for the next 1000
        mutated sequences. implemented for increased efficiency'''
        return iter(np.random.binomial(self.sequence.tag_len,self.settings['mut_rate'],1000))


    def mutate_seq(self,pop,seq_id_new,seq_id_old):
        #do the original mutation as before
        seq_sim.Simulation.mutate_seq(self,pop, seq_id_new,seq_id_old)
        tag = self.current_gen.tags[seq_id_old]
        #add new functionality: mutate tag
        try:
            nr_mut = self.mutations_per_tag.next()
        except StopIteration:
            self.mutations_per_tag = self.new_mutations_per_tag()
            nr_mut = self.mutations_per_tag.next()
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
            pos_tags = np.random.choice(range(4**self.sim.sequence.tag_len),len(self.sim.sequence.tag_dist))
            self.tags = np.random.choice(pos_tags,self.n_seq,p=self.sim.sequence.tag_dist)
        else:
            self.tags = tags

    def get_tag(self,seqID):
        return self.sim.sequence.get_tag(self.tags[seqID])

    def add_sequence(self,tag=-1,changes=None):
        new_seq_id = seq_sim.Population.add_sequence(self,changes)
        self.tags.append(tag)
        return new_seq_id



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

if __name__ == '__main__':
    tag_dist = scats.expon.pdf(range(1000),scale=400)
    np.random.shuffle(tag_dist)
    tag_dist = tag_dist/sum(tag_dist)

    sim = Simulation(tag_dist,n_seq_init=10000)
    for i in range(10):
        print '#', i
        sim.new_generation(new_gen=Population(sim,tags=[],n_seq=0))
        counts = Counter(sim.current_gen.tags)
        for i in counts:
            print i, counts[i]
