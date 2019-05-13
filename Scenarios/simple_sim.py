import SeqSimEvo as seq_sim
import argparse

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

if __name__ == '__main__':
    #parse command line arguments
    parser = argparse.ArgumentParser(description='simple simulation of sequence evolution')
    parser.add_argument('-n', nargs=1, type=int,
                        help='number of generations to run the simulation for')
    parser.add_argument('-o', nargs=1,default = 'HIV',
                        help='organism simulation settings to use')
    parser.add_argument('-st', nargs='+',default = -1,type=float,
                        help='list of sampling times, -1 for last, -0.5 for midpoint,\
                              defaults to last generation')
    parser.add_argument('-ss',nargs='+',default=10,type=int,
                        help='list of sampling sizes for each sampling point, \
                              or single value (use same sample size for all timepoints),\
                              defaults to 10 per time point')

    args = parser.parse_args()

    settings = {}
    settings['n_gen'] = args.n[0]
    settings['sample_sizes'] = args.ss


    if type(args.st) is not list:
        if args.st == -1:
            settings['sample_times'] = [args.n]
        elif args.st == -0.5:
            settings['sample_times'] = [args.n/2]
        else:
            settings['sample_times'] = [args.st]
    else:
        settings['sample_times'] =  [int(i) if i > 0 else int(args.n[0]*-i) for i in args.st ]


    #run scenario
    print run(settings, args.o)
