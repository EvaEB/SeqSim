import seq_sim
import numpy as np


class multiple_compartments(object):
    '''
    simulation of multiple compartments

    keyword arguments

    * sim_settings: simulation settings to use, as in the base 'simulation' class (default: 'HIV')

    * n_comparments: number of compartments to simulate (default: 2)

    * migration: either a float (equal migration rates between compartments), or
      a migration matrix with rows the 'from' compartment and columns the 'to'
      compartment (default 0)

    * diverse_index: differences in fitness table between the compartments. 0 is
      no differences (parallel evolution), 1 is all effects are different. (default 0)

    * names: names of the compartments (default: 0,1,2,3,...)
    '''
    def __init__(self,sim_settings = 'HIV',n_comparments = 2,migration=0.0,
                 diverse_index = 0, names=None,n_seq_init=1,**kwargs):
        self.base_sim = seq_sim.Simulation(simulation_settings=sim_settings)
        self.n_comparments = n_comparments

        #initialize migration matrix
        if type(migration) is float:
            self.mig_matrix = np.ones((n_comparments,n_comparments))*migration
        else:
            migration = np.array(migration)
            if np.shape(migration) == (n_comparments,n_comparments):
                self.mig_matrix = migration
            else:
                except_string = 'wrong dimensions of migration matrix. is {}, should be {} or float '.format(np.shape(migration),(n_comparments,n_comparments))
                raise BaseException(except_string)

        #set up sims for all compartments
        self.sims = []
        for i in range(n_comparments):
            if type(n_seq_init) is int :
                n_seq = n_seq_init
            else:
                n_seq = n_seq_init[i]


            if names is not None:
                name = names[i]
            else:
                name = str(i)

            new_kwargs = {j: kwargs[j][i] for j in kwargs}
            self.sims.append(self.base_sim.copy(name,n_seq=n_seq,**new_kwargs))


    def new_generation(self):
        #do migration
        #get emigrants from each compartment
        migrants = [[] for i in range(self.n_comparments)]
        for comp1 in range(self.n_comparments):
            for comp2 in range(self.n_comparments):
                if comp1 != comp2:
                    #select migrants
                    mig_rate = self.mig_matrix[comp1,comp2]
                    n_to_migrate = mig_rate*self.sims[comp1].current_gen.n_seq
                    if n_to_migrate < 1:
                        if np.random.random() < n_to_migrate:
                            n_to_migrate = 1
                        else:
                            n_to_migrate = 0
                    else:
                        n_to_migrate = int(n_to_migrate)
                    to_migrate = self.sims[comp1].current_gen.get_sample(n_to_migrate)

                    # remember changes in these migrants for later adding, remove from the current compartment
                    for seq in to_migrate:
                        migrants[comp2].append(self.sims[comp1].current_gen.get_seq(seq))
                        self.sims[comp1].current_gen.delete_sequence(seq)


        #add immigrants to each compartment
        for comp in range(self.n_comparments):
            for change in migrants[comp]:
                self.sims[comp].current_gen.add_sequence(changes = change)

            #let this compartment generate a new generation
            self.sims[comp].new_generation(dieout=True)
        if sum([i.n_seq for i in self.sims]) == 0:
            dieout = True
            for i in self.sims:
                i.n_seq = i.settings['n_seq_init']
            self.new_generation()


    def __str__(self):
        out = ''
        for i in self.sims:
            out+='compartment {}\t {} seqs\n'.format(i.settings['name'],i.current_gen.n_seq)
        return out

##TODO: implement temporal sampling
def run(scenario_settings,organism_settings):
    sim = multiple_compartments(sim_settings = organism_settings,
                                n_comparments = scenario_settings['n_comparments'],
                                diverse_index = scenario_settings['diverse_index'],
                                names=scenario_settings['names'],
                                n_seq_init=scenario_settings['n_seq_init'],
                                migration=scenario_settings['migration'],
                                mut_rate=scenario_settings['mut_rate'],
                                R0=scenario_settings['R0'],
                                max_pop=scenario_settings['max_pop'])

    if type(scenario_settings['sampling_amount']) is int:
        n_seq = [scenario_settings['sampling_amount']]*scenario_settings['n_comparments']
    else:
        n_seq = scenario_settings['sampling_amount']

    fasta = ''

    for i in range(scenario_settings['n_gen']):
        print i
        sim.new_generation()
        if i+1 in scenario_settings['sampling_times']:
            for j,s in enumerate(sim.sims):
                fasta += s.current_gen.to_fasta(n_seq=n_seq[j],
                                                description='-{}-gen{}'.format(s.settings['name'],i+1))

    return fasta

if __name__ == '__main__':
    sim = multiple_compartments(migration=[[0,0.005],[0.005,0]],n_seq_init=[1,0],
                                mut_rate=[0,0.0001],R0=[1.1,6],names=['blood', 'liver'],
                                max_pop=[1000,10000])
    #test max_pop
    for i in range(20):
        sim.new_generation()
        print sim

    print sim
