import seq_sim
import numpy as np
from collections import Counter
class Simulation(seq_sim.Simulation):
    def __init__(self,interactionTable=None,interactionFraction=0.3,**kwargs):
        seq_sim.Simulation.__init__(self,**kwargs)
        if interactionTable is None:
            self.interactionTable = self.__generateInteractions(fraction=interactionFraction)
        else:
            self.interactionTable = interactionTable


    def transform_fitnessTable(self,fraction=0.3,transform=['inc',0.05]):
        n_to_change = int(len(self.sequence)*3*fraction)
        #print self.sequence.sequence
        options = [(i,x) for i in range(len(self.sequence)) for x in range(4) if x != self.sequence[i]]
        indx = np.random.choice(range(len(options)),n_to_change,replace=False)
        for i in indx:
            old_value = self.fitness_table[options[i][1],options[i][0]]
            if transform[0] == 'inc':
                self.fitness_table[options[i][1],options[i][0]] = old_value + transform[1]
            elif transform[0] == 'dec':
                self.fitness_table[options[i][1],options[i][0]] = old_value - transform[1]
                if self.fitness_table[options[i][1],options[i][0]] < 0:
                    self.fitness_table[options[i][1],options[i][0]] = 0
            elif transform[0] == 'set':
                self.fitness_table[options[i][1],options[i][0]] = transform[1]
            else:
                raise ValueError('transform[0] must be one of "inc","dec","set"')

    def __generateInteractions(self,fraction,ratio=0.5):
        n_interactions = int(len(self.sequence)*(len(self.sequence)-1)*fraction)
        combis = []
        for i in range(n_interactions):
            while True:
                combi = [[np.random.randint(len(self.sequence)),np.random.randint(4)] for i in range(2)]
                if combi[0][0] != combi[1][0]:
                    break
            if combi[0][0]>combi[1][0]:
                combi = [combi[1],combi[0]]
            if np.random.random() < ratio:
                combi.append(1)
            else:
                combi.append(-1)
            combis.append(combi)
        return combis

    def get_nr_offspring(self, sequence_id, return_fitness=False):
        """returns the number of offspring of a sequence according to the fitness
        of that sequence"""
        R0 = self.settings['R0']
        changes = self.current_gen.get_seq(sequence_id)
        fitness = 1
        if changes is not None:
            for pos, base in zip(changes[:, 0], changes[:, 1]):
                fitness *= (self.fitness_table[int(base), int(pos)])
            if len(changes)>1:
                sort = np.argsort(changes[:,0])
                for i,change1 in enumerate(changes[sort]):
                    in_list = [all(i[0] == change1) for i in self.interactionTable]
                    if any(in_list):
                        interactions = [self.interactionTable[i] for i in np.where(in_list)[0]]
                        for change2 in changes[sort[1+i:]]:
                            print [all(i[1] == change2) for i in interactions]

        if return_fitness:
            return np.random.poisson(fitness*R0), fitness
        return np.random.poisson(fitness*R0)

    def copy(self,name):
        return Simulation(name=name,interactionTable=self.interactionTable)


class Population(seq_sim.Population):
    def __init__(self,**kwargs):
        seq_sim.Population.__init__(self,**kwargs)

class Seq(seq_sim.Seq):
    def __init__(self,**kwargs):
        seq_sim.Seq.__init__(self,**kwargs)



if __name__ == '__main__':
    env1 = Simulation(name='env1',interactionFraction=0.30)
#    env2 = env1.copy('env2')
#    env2.transform_fitnessTable(transform=['set',0])
    for i in range(10):
        print i
        env1.new_generation()
