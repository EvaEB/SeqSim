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

def skyline(plot=False,plot_freq=1,progress=False):
    events =     [(10,'t',0.1),(20,'v',1e5)] #(time,eventtype, new value), events: t: tranfer prop change, v: total volume change
    #event_times = [i[0] for i in events]
    initial_size = 1e4
    transfer_props = 0.01
    max_passage = 100
    skyline = passaging('phix174', initial_size,initial_size, transfer_props,2)

    if plot:
        pp = progress_plot()
    if progress:
        prog = tqdm.tqdm(range(max_passage))

    total_time = 0
    for i in tqdm.tqdm(range(max_passage)):
        while len(events) > 0 and events[0][0] == i:
            if events[0][1] == 't':
                skyline.transfer_prop = events[0][2]
            elif events[0][1] == 'v':
                skyline.max_size = events[0][2]
            del events[0]
            #del event_times[0]
        skyline.passage(1)
        total_time+=1

        if plot and (i%plot_freq==0 or i == max_passage-1):
            pp.update_plot(total_time, skyline)

    plt.ioff()
    plt.show()
    print skyline.current_gen.to_fasta(n_seq=30)

skyline(plot=True, plot_freq=5,progress=True)
