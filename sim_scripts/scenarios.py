from virus_passaging import passaging

import matplotlib.pyplot as plt
import numpy
import time
import tqdm


class progress_plot():
    def __init__(self):
        plt.ion()
        self.fig = plt.figure()
        self.ax1 = self.fig.add_subplot(2,2,1)
        self.ax2 = self.fig.add_subplot(2,2,2)
        self.ax3 = self.fig.add_subplot(2,2,3)
        self.ax4 = self.fig.add_subplot(2,2,4)

        self.hl, = self.ax1.plot([],[])
        self.ax1.set_title('population size')

        self.frac, = self.ax2.plot([],[])
        self.ax2.set_title('maximum fraction of a mutation')

        self.mut, = self.ax3.plot([],[])
        self.ax3.set_title('total mutations')

        self.unique, = self.ax4.plot([],[])
        self.ax4.set_title('unique mutations')

    def update_line(self,fig,line, x,y):
        line.set_xdata(numpy.append(line.get_xdata(), x))
        line.set_ydata(numpy.append(line.get_ydata(), y))
        fig.canvas.draw()
        fig.canvas.flush_events()

    def update_plot(self,time,sim):
        self.update_line(self.fig,self.hl,time,sim.current_gen.n_seq)

        stats = sim.current_gen.stats()
        self.update_line(self.fig,self.frac,time,stats['max_fraction'])
        self.update_line(self.fig,self.mut,time,stats['total_mutations'])
        self.update_line(self.fig,self.unique,time,stats['unique_mutations'])

        self.ax1.relim()
        self.ax1.autoscale_view()
        self.ax2.relim()
        self.ax2.autoscale_view()
        self.ax3.relim()
        self.ax3.autoscale_view()
        self.ax4.relim()
        self.ax4.autoscale_view()


def run(scenario,scenario_settings,organism_settings):
    scenario = scenario.split('_')[0]
    if scenario == 'SimpleSim':
        import simple_sim
        fasta = simple_sim.run(scenario_settings,organism_settings)

    elif scenario == 'VirusPassaging':
        import virus_passaging
        fasta = virus_passaging.run(scenario_settings,organism_settings)
    elif scenario == 'RecreateDataset':
        import recreate_dataset
        fasta = recreate_dataset.recreate_dataset(scenario_settings['SampleSize'],
                                                  scenario_settings['nrMutations'],
                                                  scenario_settings['apobec_IDs'],
                                                  apobec_rate=scenario_settings['apobec_rate'],
                                                  action='fasta',
                                                  simulation_settings=organism_settings)
    elif scenario == 'MultipleCompartments':
        import multiple_compartments
        fasta = multiple_compartments.run(scenario_settings,organism_settings)

    return fasta

if __name__ == '__main__':
    #events =     [(10,'t',0.1),(20,'v',1e5)] #(time,eventtype, new value), events: t: tranfer prop change, v: total volume change
    #skyline(events, plot=True, plot_freq=5,progress=True)
    control(plot=True,plot_freq=10)
    #migration()
    #root()
