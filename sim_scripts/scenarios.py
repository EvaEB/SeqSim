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

def skyline(events,initial_size=1e4,transfer_props=0.01,max_passage=100,
            gen_per_transfer=2,plot=False,plot_freq=1,progress=False):
    skyline_sim = passaging('phix174', initial_size,initial_size, transfer_props,
                            gen_per_transfer)
    print skyline_sim.current_gen.to_fasta(seq_ids=[0],description='consensus'),
    if plot:
        pp = progress_plot()

    total_time = 0
    for i in tqdm.tqdm(range(max_passage)):
        while len(events) > 0 and events[0][0] == i:
            if events[0][1] == 't':
                skyline_sim.transfer_prop = events[0][2]
            elif events[0][1] == 'v':
                skyline_sim.max_size = events[0][2]
            del events[0]
            #del event_times[0]
        skyline_sim.passage(1)
        total_time+=1

        if plot and (i%plot_freq==0 or i == max_passage-1):
            pp.update_plot(total_time, skyline_sim)

    plt.ioff()
    plt.show()
    print skyline_sim.current_gen.to_fasta(n_seq=30)

def control(initial_size=1e4,transfer_props=0.1,max_passage=100,
            gen_per_transfer=2,plot=False,plot_freq=1,progress=False):
    control_sim = passaging('phix174', initial_size,initial_size, transfer_props,
                            gen_per_transfer)
    print control_sim.current_gen.to_fasta(seq_ids=[0],description='consensus')
    if plot:
        pp = progress_plot()

    total_time = 0
    for i in tqdm.tqdm(range(max_passage)):
        control_sim.passage(1)
        total_time+=1

        if plot and (i%plot_freq==0 or i == max_passage-1):
            pp.update_plot(total_time, control_sim)

    plt.ioff()
    plt.show()
    print control_sim.current_gen.to_fasta(n_seq=30)

def migration():
    transfer_props = [[1e-2, 0, 5e-3],
                      [5e-3, 1e-2, 0],
                      [5e-3, 0, 1e-2]]
    max_passage =  100
    n_pop = len(transfer_props)
    n_seq_init = 1e4
    viral_pops = [passaging('phix174',initial_size=n_seq_init,
                            transfer_prop = transfer_props[0][0])]
    print viral_pops[0].current_gen.to_fasta(seq_ids=[0],description='consensus')

    for i in range(1,n_pop):
        viral_pops.append(viral_pops[-1].copy(i,n_seq=n_seq_init))
        viral_pops[-1].transfer_prop = transfer_props[i][i]

    for i in tqdm.tqdm(range(max_passage)):
        for i in viral_pops:
            i.passage(1)
        for i in range(n_pop):
            for j in range(n_pop):
                if i != j:
                    amount = transfer_props[i][j]
                    migration_size = int(viral_pops[j].current_gen.n_seq*amount)
                    migrating = viral_pops[j].current_gen.get_sample(migration_size)
                    for m in migrating:
                        changes = viral_pops[j].current_gen.get_seq(m)
                        viral_pops[i].current_gen.add_sequence(changes)

    for i,pop in enumerate(viral_pops):
        print pop.current_gen.to_fasta(n_seq=10,description='_population_'+str(i))

def root():
    times = [[10,0],[20,1],[25,2],[30,3]]
    total_time = 35
    times.append([total_time+1,1])
    n_seq_init = 100
    viral_pops = viral_pops = [passaging('phix174',initial_size=n_seq_init)]

    print viral_pops[0].current_gen.to_fasta(seq_ids=[0],description='consensus')

    for i in tqdm.tqdm(range(total_time)):
        if i == times[0][0]:
            viral_pops.append(viral_pops[times[0][1]].copy(i,n_seq=n_seq_init))
            times.pop(0)
        for j in viral_pops:
            j.passage(1)

    for i,pop in enumerate(viral_pops):
        print pop.current_gen.to_fasta(n_seq=10,description='_population_'+str(i))



if __name__ == '__main__':
    #events =     [(10,'t',0.1),(20,'v',1e5)] #(time,eventtype, new value), events: t: tranfer prop change, v: total volume change
    #skyline(events, plot=True, plot_freq=5,progress=True)
    #control(plot=True,plot_freq=10)
    #migration()
    root()
