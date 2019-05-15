import SeqSimEvo as seq_sim
import numpy as np
from collections import Counter
from copy import deepcopy
from tqdm import tqdm

class Simulation(seq_sim.Simulation):
    def __init__(self,interactionTable=None,interactionFraction=0.3, ratio=0.5,**kwargs):
        seq_sim.Simulation.__init__(self,**kwargs)
        if interactionTable is None:
            self.interactionTable = self.generateInteractions(fraction=interactionFraction, ratio=ratio)
        else:
            self.interactionTable = interactionTable


    def transform_fitnessTable(self,fraction=0.3,transform=['inc',0.05]):
        n_to_change = int(len(self.sequence)*3*fraction)
        #print self.sequence.sequence
        options = [(i,x) for i in range(len(self.sequence)) for x in range(4) if x != self.sequence[i]]
        indx = np.random.choice(range(len(options)),n_to_change,replace=False)
        for i in indx:
            old_value = self.fitness_table[options[i][1],options[i][0]]
            if transform[0] == 'product':
                self.fitness_table[options[i][1],options[i][0]] = old_value * transform[1]
            elif transform[0] == 'add':
                self.fitness_table[options[i][1],options[i][0]] = old_value + transform[1]
                if self.fitness_table[options[i][1],options[i][0]] < 0:
                    self.fitness_table[options[i][1],options[i][0]] = 0
            elif transform[0] == 'set':
                self.fitness_table[options[i][1],options[i][0]] = transform[1]
            elif transform[0] == 'scramble':
                self.fitness_table[options[i][1],options[i][0]] = np.random.choice(self.fitness_table.flatten())
            else:
                raise ValueError('transform[0] must be one of "product","add","set","scramble"')

    def generateInteractions(self,fraction,ratio=0.5):
        n_interactions = int(len(self.sequence)*(len(self.sequence)-1)*fraction)
        combis = {}
        for i in tqdm(range(n_interactions)):
            while True:
                combi = [[np.random.randint(len(self.sequence)),np.random.randint(4)] for i in range(2)]
                if combi[0][0] != combi[1][0]:
                    break
            if combi[0][0]>combi[1][0]:
                combi = [combi[1],combi[0]]
            if np.random.random() < ratio:
                #combi.append(1)
                fit1 = self.get_fitness_effect(combi[0][0],combi[0][1])
                fit2 = self.get_fitness_effect(combi[1][0],combi[1][1])
                cur = fit1*fit2
                new = 1+np.random.exponential(self.settings['parameters']['lb'])
                if new>cur:
                    combi.append(new)
                else:
                    combi.append(cur)
            else:
                combi.append(0)
            try:
                combis[tuple(combi[0])].append(combi[1:])
            except KeyError:
                combis[tuple(combi[0])] = [combi[1:]]
        return combis

    def get_nr_offspring(self, sequence_id, return_fitness=False):
        """returns the number of offspring of a sequence according to the fitness
        of that sequence"""
        R0 = self.settings['R0']
        changes = self.current_gen.get_seq(sequence_id)
        fitness = 1
        if changes is not None:
            #for pos, base in zip(changes[:, 0], changes[:, 1]):
            #    fitness *= (self.fitness_table[int(base), int(pos)])
            sort = np.argsort(changes[:,0])
            already_done = []
            for i,change1 in enumerate(changes[sort]):
                if tuple(change1) in self.interactionTable:
                    ch1 = [tuple(c[0]) for c in self.interactionTable[tuple(change1)]]
                    chs2 = set([tuple(c) for c in changes[sort][i+1:]])
                    interactions = set(ch1).intersection(chs2)
                    for interaction in interactions:
                        already_done += [tuple(change1),interaction]
                        value = [self.interactionTable[tuple(change1)][i][-1] for i in range(len(ch1)) if ch1[i]==interaction]
                        fitness*=np.product(value)
                done=False
                if len(already_done)>0:
                    if tuple(change1) in already_done:
                        done=True
                if not done:
                    fitness*=(self.fitness_table[int(change1[1]), int(change1[0])])

        if return_fitness:
            return np.random.poisson(fitness*R0), fitness
        return np.random.poisson(fitness*R0)

    def copy(self,name):
        return Simulation(simulation_settings=deepcopy(self.settings), sequence = self.sequence,
                          fitness_table=deepcopy(self.fitness_table),name= name,
                          interactionTable=deepcopy(self.interactionTable))


class Population(seq_sim.Population):
    def __init__(self,**kwargs):
        seq_sim.Population.__init__(self,**kwargs)

class Seq(seq_sim.Seq):
    def __init__(self,**kwargs):
        seq_sim.Seq.__init__(self,**kwargs)

def analyze_passage(sim,f,passage):
    try:
        changes = np.vstack(list(sim.current_gen.changes.values()))
    except ValueError:
        return 0
    n_seq = float(sim.n_seq)

    changes = [tuple(i) for i in changes]
    counts = Counter(changes)
    for i in counts:
        f.write("{}\t{}\t{}\t{}\t{}\n".format(i[0], i[1], counts[i]/n_seq,
                                              passage,sim.settings['name']))

def run_custom():
    interactionfile = sys.argv[1]
    resultfile = sys.argv[2]
    timingsfile = sys.argv[3]

    interactionFraction = float(sys.argv[4])
    ratio = float(sys.argv[5])
    transformType = sys.argv[6]
    transformValue = float(sys.argv[7])
    fraction = float(sys.argv[8])

    env1_1 = Simulation(name='env1_1',interactionFraction=interactionFraction,
                        ratio=ratio,seq_len=500)
    sims = [env1_1]+[env1_1.copy('env1_'+str(i)) for i in range(2,6)]
    env2_1 = env1_1.copy('env2_1')
    env2_1.transform_fitnessTable(transform=[transformType,transformValue],fraction=fraction)

    sims += [env2_1]+[env2_1.copy('env2_'+str(i)) for i in range(2,6)]

    with open(interactionfile,'w') as f:
        pickle.dump(env1_1.interactionTable,f)

    f = open(resultfile,'w',buffering=1)
    t = open(timingsfile,'w',buffering=1)

    for i in range(201):
        for sim in sims:
            if i%2 == 0:
                sim.settings['max_pop'] = 100
            sim.new_generation()
            sim.settings['max_pop'] = 100000

        if i%10 == 0:
            t.write("{} {}\n".format(time.time(), i))
            for sim in sims:
                analyze_passage(sim, f,i)
    f.close()
    t.close()

if __name__ == '__main__':
    import time
    import pickle
    import sys

    interactionfile = sys.argv[1]
    resultfile = sys.argv[2]
    timingsfile = sys.argv[3]

    env1_1 = Simulation(name='env1_1',interactionFraction=0.05,
                        ratio=0,seq_len=500)
    sims = [env1_1]+[env1_1.copy('env1_'+str(i)) for i in range(2,6)]
    env2_1 = env1_1.copy('env2_1')
    #nv2_1.generateInteractions(float(sys.argv[4]),1)
    #env2_1.transform_fitnessTable(transform=['scramble',0],fraction=0.3)
    #env2_1.transform_fitnessTable(transform=['product',1.1],fraction=1)
    #env2_1.transform_fitnessTable(transform=['set',0],fraction=0.3)
    #env2_1.transform_fitnessTable(transform=['set',1],fraction=0.3)



    sims += [env2_1]+[env2_1.copy('env2_'+str(i)) for i in range(2,6)]

    with open(interactionfile,'wb') as f:
        pickle.dump(env1_1.interactionTable,f)

    f = open(resultfile,'w',buffering=1)
    t = open(timingsfile,'w',buffering=1)

    for i in range(10):
        for sim in sims:
            if i%2 == 0:
                if 'env2' in sim.settings['name']:
                    sim.settings['max_pop'] = int(sys.argv[4])
                else:
                    sim.settings['max_pop'] = 1000

            sim.new_generation()
            sim.settings['max_pop'] = 10000

        if i%10 == 0:
            t.write("{} {}\n".format(time.time(), i))
            for sim in sims:
                analyze_passage(sim, f,i)
    f.close()
    t.close()
