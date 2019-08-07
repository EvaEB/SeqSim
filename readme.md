(c) 2019 Eva Bons

# seq_sim - simulations of viral diversification

## installation
First, make sure you have a working python installation (2.7/3.6), with the following packages installed
* numpy
* scipy
* pyyaml
* Tqdm
* Ete3
* Ipython
* matplotlib

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

# Scenarios - simulation of specific scenarios using SeqSimEvo
## The gui
### Installation
Install dependencies:
  * appJar `pip install appjar`
  * pyyaml `pip install pyyaml`
  * fasta_tools (available on my git - ask for access if you don't have it)

make sure seq_sim is in the same folder as sim_scipts (otherwise some file references don't work)

### Usage
* start `gui.py` (`python gui.py` in a terminal)
* choose the organism and scenarios
* change settings if required (changes are not saved to the original file!)
* run a simulation
* view results (fasta, highlighter or neighbor-joining tree)

### Adding new scenarios
Write a python function for the new scenario, make sure it can take a library with scenario_settings and organism_settings as arguments and returns a string formatted as a fasta-file

    new_scenario(scenario_settings, organism_settings)

Add an entry to the `run` funtion in `scenarios.py`

    if scenario == 'NewScenario':
      fasta = new_scenario(scenario_settings, organism_settings)

Add a settings file (or template) to settings_files. Filetype must be of the format `NewScenario[_something]`


### Making new settings files
 * Everything before the underscore must match what is checked in the if-statement in the `run`-function in `scenarios.py` for this scenarios
 * The file must be in the 'yaml' format (http://www.yaml.org/spec/1.2/spec.html)
 * Due to issues with scientific notation parsing, it is important to always use a decimal point when using scientific notation of numbers (`1.0e5` rather than `1e5`)

## recreate dataset
run a simulation of sequence diversification until the number of unique mutations matches that in a given dataset
*more info on command-line usage coming*

## virus passaging
run a simulation of virus passaging experiments
*more info on command-line usage coming*
