import seq_sim as sim
import numpy as np

master_sim = sim.Simulation(simulation_settings='HIV')
n_gen = 100

#settings of the different compartments
compartment_names = ['sanctuary', 'liver','blood']
compartment_sizes = [1e4,         1e6,     1e10]
initial_sizes =     [0,           0,       10]
migration_rates =  np.array([[0,           0,       1e-3],   #row=from, col = to
                             [0,           0,       1e-2],
                             [1e-3,        1e-2,    0]])
replication     =   [0.5,         1,       0.001]


#setup different compartments:
compartments = []
for i in range(len(compartment_names)):
    compartments.append(master_sim.copy(compartment_names[i], n_seq=initial_sizes[i]))
    compartments[-1].settings['max_pop'] = compartment_sizes[i]
    compartments[-1].settings['mut_rate'] = master_sim.settings['mut_rate']*replication[i]
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
    if gen == 50:
        compartments[1].settings['max_pop'] = 0

    if gen == 51:
        compartments[1].settings['max_pop'] = compartment_sizes[1]

    #print some stats
    print 'gen {}: nseq'.format(gen),
    for comp in compartments:
        print comp.current_gen.stats()['n_seq'],
    print ''
