import seq_sim

def run(scenario_settings, sim_settings):
    sim = seq_sim.Simulation(simulation_settings=sim_settings)

    if type(scenario_settings['sample_sizes']) is int:
        n_seq = [scenario_settings['sample_sizes']]*len(scenario_settings['sample_times'])
    else:
        n_seq = scenario_settings['sample_sizes']

    output = ''
    for i in range(scenario_settings['n_gen']):
        sim.new_generation()
        if i+1 in scenario_settings['sample_times']:
            index = [j for j in range(len(scenario_settings['sample_times'])) if
                        scenario_settings['sample_times'][j] == i+1][0]
            output += sim.current_gen.to_fasta(n_seq=n_seq[index],
                                               description='-gen{}'.format(i+1))

    return output
