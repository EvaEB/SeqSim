# seq_sim - simulations of viral diversification

## installation
First, make sure you have a working python installation (2.7), with the following packages installed
* numpy
* pyyaml
* progressbar

Currently, I have not managed to find an easy way to install/manage packages, sorry!
To install, you'll  have to manually put the `seq_sim` folder somewhere in your `PYTHONPATH`

## usage
### initialize the simulation
    import seq_sim

    sim = seq_sim.Simulation(simulation_settings)

simulation settings can be set by
* providing a default (currently, 'HIV' and 'phix174' are avaiable)

		sim = seq_sim.Simulation(simulation_settings)
* pointing to your own settings file (a template is available in the folder `simulation_settings`):

		sim = seq_sim.Simulation(simulation_settings = 'path/to/my/settings')
* providing the settings a keyword arguments (this is always possible, and will override settings provided in a settings file): e.g.

	 	sim = seq_sim.Simulation(mut_rate = 1e-10, model='lognormal', parameters = {'fl': 0.2, 'mu':0.2, 'sigma': 0.2})
	 	sim = seq_sim.Simulation('phix174', model='spikes')

see `help(seq_sim.Simulation)` for all possible settings

### run a simulation
    for gen in range(nr_gen):
	    sim.new_generation()

### get output
* for fasta output of any number of sequences (use `nseq = sim.current_gen.n_seq` for all sequences in the simulation)

		sim.current_gen.to_fasta(n_seq)

* to see some statistics

		sim.current_gen.stats()

# sim_scripts - specific simulation scenarios using seq_sim
## recreate dataset


## virus passaging

## using the gui
