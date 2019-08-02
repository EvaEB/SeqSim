'''
Simulation extension to implement reproduction via binary fission

In every generation, cells die or reproduce. R0 indicated the relative chance to
die or reproduce: an R0 of 1 equals initial equal chances to die or reproduce,
R0 of 2 is twice as likely to reproduce as to die. fitness values modify this
probability, as f*R0

All other functionalities (calculation of fitness, mutating sequences, etc.) are
as in the original simulation code

'''
import simulation as sim
from simulation import Population, Seq
import random

class Simulation(sim.Simulation):
    def survive(self, fitness):
        '''decides whether a sequence with a certain fitness survives or not'''
        p = (1-(1.0/(1+fitness*self.settings['R0']))) #probability of survival = f*R0/(1+f*R0)
        if random.random()<p:
            return True
        else:
            return False

    def new_generation(self):
        #Updates current_gen, effective_pop, gen, average_fitness,
        #and n_seq
        self.effective_pop = 0
        self.gen+=1

        to_delete = []
        all_fitnesses = []
        for cur_seq in range(len(self.current_gen)):
            _, fitness = self.get_nr_offspring(cur_seq,return_fitness=True)
            all_fitnesses.append(fitness)

            if self.survive(fitness):
                changes = self.current_gen.get_seq(cur_seq)
                new_index = self.current_gen.add_sequence(changes)
                mutated = self.mutate_seq(self.current_gen,cur_seq,new_index)
                self.effective_pop += (not changes is None) + mutated
            else:
                to_delete.append(cur_seq)

        for seq in sorted(to_delete,reverse=True): #delete last first so that new seqID problems arise
            self.current_gen.delete_sequence(seq)

        n_seq = len(self.current_gen)
        if n_seq > self.settings['max_pop']:
            for delete in range(n_seq-self.settings['max_pop']):
                to_del = random.randint(0,n_seq-delete-1)
                self.current_gen.delete_sequence(to_del)

        elif n_seq == 0:
            print 'died out'
            raise NotImplementedError('population died out')

        self.n_seq = len(self.current_gen)
        self.average_fitness = sum(all_fitnesses)/len(all_fitnesses)

if __name__ == '__main__':
    from tqdm import tqdm

    sim_test = Simulation(R0 = 2,max_pop=10000,n_seq_init=1000)
    for i in tqdm(range(50)):
        sim_test.new_generation()
    print sim_test.current_gen
