import SeqSimEvo as sim
import numpy as np

master_sim = sim.Simulation(simulation_settings='HIV')
n_gen = 10

#settings of the different compartments
compartment_names = ['sanctuary', 'liver','blood']
compartment_sizes = [1e3,         1e5,     1e2]
initial_sizes =     [0,           0,       10]
migration_rates =  np.array([[0,           0,       1e-3],   #row=from, col = to
                             [0,           0,       1e-2],
                             [1e-3,        1e-2,    0]])
replication     =   [0.2,         1,       0.166]
mut_rate =          [0.2,         1,       0]

fasta = ''
#setup different compartments:
compartments = []
for i in range(len(compartment_names)):
    compartments.append(master_sim.copy(compartment_names[i], n_seq=initial_sizes[i]))
    compartments[-1].settings['max_pop'] = compartment_sizes[i]
    compartments[-1].settings['mut_rate'] = master_sim.settings['mut_rate']*mut_rate[i]
    compartments[-1].settings['R0'] = master_sim.settings['R0']*replication[i]


#do simulations
for gen in range(n_gen):
    #new generation
    for comp in compartments:
        comp.new_generation(dieout=True)

    #do migration
    additions = [[] for i in compartment_names]
    for comp1 in range(len(compartments)):
        for comp2 in range(len(compartments)):
            migration = migration_rates[comp1, comp2]
            if migration>0:
                n_migrate = np.random.poisson(migration*len(compartments[comp1].current_gen))
                for i in range(n_migrate):
                    migrant = np.random.randint(len(compartments[comp1].current_gen))
                    additions[comp2].append(compartments[comp1].current_gen.get_seq(migrant))
                    compartments[comp1].current_gen.delete_sequence(migrant)
    for i, add in enumerate(additions):
        for seq in add:
            compartments[i].current_gen.add_sequence(changes = seq)

    #liver transplant
    # if gen == 50:
    #     compartments[1].settings['max_pop'] = 0
    #
    # if gen == 51:
    #     compartments[1].settings['max_pop'] = compartment_sizes[1]

    #sampling
    if gen in [4,20,50,99]:#55,70,99]:
        fasta += compartments[0].current_gen.to_fasta(n_seq=10,description='-1-{}'.format(gen))
        fasta += compartments[1].current_gen.to_fasta(n_seq=10,description='-2-{}'.format(gen))
        #fasta += compartments[2].current_gen.to_fasta(n_seq=10,description='-3-{}'.format(gen))

    #print some stats
    print 'gen {}: nseq'.format(gen),
    for comp in compartments:
        print comp.current_gen.stats()['n_seq'],
    print ''

with open('3comp.fasta', 'w') as f:
    f.write(fasta)
