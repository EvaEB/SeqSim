(c) 2019 Eva Bons

# SeqSimEvo - simulations of sequence evolution with selection

## Installation
installation requires a working installation of python (v2.7 or higher)

    python setup.py install

## Usage
Simulations can be run in several ways:
  * run a pre-existing scenario

### Initialize the simulation
    from SeqSimEvo import simulation

    sim = simulation.Simulation(simulation_settings)

simulation settings can be set by
* providing a default (currently, 'HIV' and 'phix174' are available)

		sim = seq_sim.Simulation(simulation_settings)
* pointing to your own settings file (a template is available in the folder `simulation_settings`):

		sim = seq_sim.Simulation(simulation_settings = 'path/to/my/settings')
* providing the settings as keyword arguments (this is always possible, and will override settings provided in a settings file): e.g.

	 	sim = seq_sim.Simulation(mut_rate = 1e-10, model='lognormal', parameters = {'fl': 0.2, 'mu':0.2, 'sigma': 0.2})
	 	sim = seq_sim.Simulation('phix174', model='spikes')

see `help(simulation.Simulation)` for all possible settings

### Run a simulation
    for gen in range(nr_gen):
	    sim.new_generation()

### Get output
* for fasta output of any number of sequences (use `nseq = sim.current_gen.n_seq` for all sequences in the simulation)

		sim.current_gen.to_fasta(n_seq)

* to see some statistics

		sim.current_gen.stats()

# Scenarios - simulation of specific scenarios using SeqSimEvo
## The GUI

### Usage
* run `SeqSimEvo_gui` in the terminal
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

Add a settings file (or template) to settings_files. the File type must be of the format `NewScenario[_something]`


### Making new settings files
 * Everything before the underscore must match what is checked in the if-statement in the `run`-function in `scenarios.py` for this scenarios
 * The file must be in the 'yaml' format (http://www.yaml.org/spec/1.2/spec.html)
 * Due to issues with scientific notation parsing, it is important to always use a decimal point when using scientific notation of numbers (`1.0e5` rather than `1e5`)

## Simple simulation
simply run a simulation for a certain number of generations.

run `SeqSimEvo_simpleSim [options]` in the terminal.

see `SeqSimEvo_simpleSim -h`

## recreate dataset
run a simulation of sequence diversification until the number of unique mutations matches that in a given dataset

run `SeqSimEvo_recreateDataset [options]` in the terminal.

see `SeqSimEvo_recreateDataset -h` for all options.

## multiple multiple_compartments
run a simulations with multiple compartments, with migration between the compartments

run `SeqSimEvo_multipleCompartments [options]` in the terminal.

see `SeqSimEvo_multipleCompartments -h`
