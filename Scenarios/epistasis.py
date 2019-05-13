import SeqSimEvo as seq_sim
import numpy as np
from collections import Counter
from copy import deepcopy
from tqdm import tqdm
import sys

class Simulation(seq_sim.Simulation):
    def __init__(self,epistasis={}, ratio=0.5,fitness_table=None,**kwargs):
        seq_sim.Simulation.__init__(self,**kwargs)
        if fitness_table is None:
            self.fitness_table = self.__getfitness_table()
        else:
            self.fitness_table = fitness_table
        self.settings['epistasis'] = epistasis

    def __getfitness_table(self):
        return {}

    def assign_value(self,changes):
        if len(changes) == 1: #single mutation: draw from MFED
            if self.settings['model'] == 'exponential':
                pars = self.settings['parameters']
                random_nr = np.random.random()
                if random_nr < pars['fl']: #lethal mutation
                    return 0
                elif random_nr < (pars['fl']+pars['fd']):
                    fit = 1-np.random.exponential(pars['ld'])
                    if fit<0:
                        fit=0
                    return fit
                elif random_nr > (pars['fl']+pars['fd']+pars['fn']):
                    return 1+np.random.exponential(pars['lb'])
            if self.settings['model'] == 'neutral':
                return 1
            if self.settings['model'] == 'lognormal':
                random_nr = np.random.random()
                pars = self.settings['parameters']
                if random_nr < pars['fl']: #lethal mutation
                    return 0
                else:
                    return np.random.lognormal(mean=pars['mu'],sigma=pars['sigma'])
            else:
                raise('model {} not implemented'.format(self.settings['model']))
        else:
            epi = self.settings['epistasis']
            random_nr = np.random.random()
            if epi['model'] == 'exponential':
                random_nr = np.random.random()
                if random_nr < epi['par']['fl']:
                    e = 0
                elif random_nr < epi['par']['fl']+epi['par']['fd']:
                    e = np.exp(-np.random.exponential(epi['par']['ld']))
                else:
                    e = np.exp(np.random.exponential(epi['par']['lb']))
            if epi['model'] == 'none':
                e = 1
            else:
                print 'model {} not implemented'.format(epi['model'])
                raise()

            if len(changes) == 2:
                str1 = '{}-{}'.format(changes[0][0],changes[0][1])
                str2 = '{}-{}'.format(changes[1][0],changes[1][1])

                try:
                    f1 = self.fitness_table[str1]
                except KeyError:
                    self.fitness_table[str1] = self.assign_value([changes[0]])
                    f1 = self.fitness_table[str1]

                try:
                    f2 = self.fitness_table[str2]
                except KeyError:
                    self.fitness_table[str2] = self.assign_value([changes[1]])
                    f2 = self.fitness_table[str2]
            else:
                found=False
                for c in range(len(changes)):
                    subset = np.vstack([changes[:c],changes[c+1:]])
                    changes_str = ['{}-{}'.format(i[0],i[1]) for i in subset]
                    changes_str = ' '.join(sorted(changes_str))
                    if changes_str in self.fitness_table.keys():
                        f1 = self.fitness_table[changes_str]
                        found=True
                        break
                if not found:
                    self.fitness_table[changes_str] = self.assign_value(subset)
                    f1 = self.fitness_table[changes_str]
                try:
                    f2 = self.fitness_table['{}-{}'.format(changes[c][0],changes[c][1])]
                except KeyError:
                    self.fitness_table['{}-{}'.format(changes[c][0],changes[c][1])] = self.assign_value([changes[c]])
                    f2 = self.fitness_table['{}-{}'.format(changes[c][0],changes[c][1])]
        return e*f1*f2

    def get_nr_offspring(self, sequence_id, return_fitness=False):
        """returns the number of offspring of a sequence according to the fitness
        of that sequence"""
        R0 = self.settings['R0']
        changes = self.current_gen.get_seq(sequence_id)
        if changes is None:
            fitness = 1
        else:
            changes_str = ['{}-{}'.format(i[0],i[1]) for i in changes]
            changes_str = ' '.join(sorted(changes_str))
            try:
                fitness = self.fitness_table[changes_str]
            except KeyError:
                self.fitness_table[changes_str] = self.assign_value(changes)
                fitness = self.fitness_table[changes_str]
        if return_fitness:
            return np.random.poisson(R0*fitness), fitness
        else:
            return np.random.poisson(R0*fitness)

    def copy(self,name):

        return Simulation(simulation_settings=deepcopy(self.settings),
                          epistasis=self.settings['epistasis'],
                          sequence = self.sequence,
                          fitness_table=deepcopy(self.fitness_table),
                          name= name)


class Population(seq_sim.Population):
    def __init__(self,**kwargs):
        seq_sim.Population.__init__(self,**kwargs)


class Seq(seq_sim.Seq):
    def __init__(self,**kwargs):
        seq_sim.Seq.__init__(self,**kwargs)


if __name__ == '__main__':
    n_gen = 200

    pars = {'fl': 0.1,'mu': -0.25,'sigma': 0.15}
    epistasis = {'same':1, 'other':0,
                 'par':{'fb':0,'lb':0,'fl':0,'fd':0,'ld':0,'fn':1},
                 'model':'none'}
    sim = Simulation(model='lognormal',
                     parameters = pars,
                     epistasis=epistasis,
                     n_seq_init=1000,
                     max_pop=10000,
                     mut_rate=1e-5,
                     seq_len=1000)

    for gen in range(n_gen):
        sim.new_generation()
        if gen%10 == 0:
            print 'mean: {}'.format(sim.average_fitness)
