import numpy as np
from copy import deepcopy
from simulation import Simulation, Seq


class Population():
    def __init__(self, simulation, n_seq=1):
                            #    id,changed,parent,number
        self.changes = np.array([[0, -1,     -1,    1    ]])

        self.seq_ids = [0]
        self.sim = simulation
        #self.sim = simulation
        self.curID = len(self.changes)-1
    def __len__(self):
        return sum(self.changes[:,3])

    def __str__(self):
        string = '#mutID (from-pos-to)\tsequence\tpatient\n'
        for i,seq in enumerate(self.seq_ids):
            parent = self.changes[seq,2]
            cur_line = seq
            changes = []
            while parent != -1:
                changes.append(self.changes[cur_line,1])
                cur_line = int(parent)
                parent = self.changes[cur_line,2]
            for change in changes:
                pos = int(change)
                new = int((change-pos)*10)
                string+='{orig}-{pos}-{to}\t{seq}\t{patient}\n'.format(orig=self.sim.sequence[pos],
                                                               pos=pos,
                                                               to=new,
                                                               seq=i,
                                                               patient=self.sim.settings['name'])
        return string

    def copy():
        pass
        ##still to implement

    def next_ID(self):
        self.curID +=1
        return self.curID


    def print_sample(self, seq_ids):
        string = '#mutID (from-pos-to)\tsequence\tpatient\n'
        for i in seq_ids:
            seq = self.seq_ids[i]
            parent = self.changes[seq,2]
            cur_line = seq
            changes = []
            while parent != -1:
                changes.append(self.changes[cur_line,1])
                cur_line = int(parent)
                parent = self.changes[cur_line,2]
            for change in changes:
                pos = int(change)
                new = int((change-pos)*10)
                string+='{orig}-{pos}-{to}\t{seq}\t{patient}\n'.format(orig=self.sim.sequence[pos],
                                                               pos=pos,
                                                               to=new,
                                                               seq=i,
                                                               patient=self.sim.settings['name'])
        print string

    def get_sample(self, sample_size):
        try:
            return np.random.choice(self.seq_ids,size=int(sample_size),replace=False)
        except ValueError:
            return self.seq_ids

    def add_sequence(self, changes=None):

        parent = 0
        if changes is not None:
            changes = sorted(changes)

            for i,change in enumerate(changes):
                parent_pos = (self.changes[:,2] == parent) & (self.changes[:,1] == change)
                if sum(parent_pos)!=0:
                    parent = np.where(parent_pos)[0][0]
                else:
                    loose_nodes = set.difference(set(np.where(self.changes[:,3] == 0)[0]), np.append(self.changes[:,2], [parent]))
                    if len(loose_nodes) == 0:
                        self.changes = np.append(self.changes,[[self.next_ID(),change,parent,0]],axis=0)
                        parent = self.curID
                    else:
                        cur_node = list(loose_nodes)[0]
                        self.changes[cur_node,:] = [cur_node,change,parent,0]
                        parent =cur_node
        self.changes[parent,3]+=1
        self.seq_ids.append(parent)
        # return self.n_seq-1

    #def add_change(self, seq_id, pos, target):

    def stats(self):
        print self.changes
        stats = {}
        stats['n_seq'] = len(self.seq_ids)
        stats['unmutated'] = self.changes[0,3]


        # stats['total_mutations'] = 0
        # for i,node in enumerate(self.changes[1:]):
        #     if node[3] > 0:
        #         n_muts = 1
        #         parent = int(node[2])
        #         while parent!=0:
        #             parent = int(self.changes[parent,2])
        #             mut_counts[parent-1]+=1
        #             n_muts+=1
        #         stats['total_mutations'] += n_muts*node[3]

        stats['unique_mutations'] = len(set(self.changes[1:,1]))

        # print mut_counts
        return stats


        #     mut_counts = np.array(Counter(all_mutations).values())
        #     if len(mut_counts) > 0:
        #         stats['majority_mutations'] = sum(mut_counts > (stats['n_seq']/2.0))
        #         stats['max_fraction'] = max(mut_counts/float(stats['n_seq']))
        #     else:
        #         stats['majority_mutations'] = 0
        #         stats['max_fraction'] = 0
        #     GA = 0
        #     for i in all_mutations:
        #         if self.sim.sequence[i[0]]==1 and i[1] == 0:
        #             GA+=1.0
        #     try:
        #         stats['GA_rate'] = GA/len(all_mutations)
        #     except ZeroDivisionError:
        #         stats['GA_rate'] = None
        #     return stats
        #
        #
        # def to_fasta(self, seq_ids=[], n_seq=None, description='',progress=False):
        #     string = ''
        #     if len(seq_ids) == 0:
        #         if n_seq is None:
        #             n_seq = self.n_seq
        #         seq_ids = random.sample(range(self.n_seq), n_seq)
        #     for i in range(len(seq_ids)):
        #         seqID = seq_ids[i]
        #         string += '>'+str(seqID)+''+str(description)+'\n'
        #         changed_here = self.get_seq(seqID)
        #         seq = deepcopy(self.sim.sequence)
        #         if changed_here is not None:
        #             for i in changed_here:
        #                 seq.sequence[int(i[0])] = int(i[1])
        #         string += str(seq)+'\n'
        #     return string
        #
        # def consensus_sequence(self):
        #     seq = deepcopy(self.sim.sequence)
        #     all_mutations = np.vstack(self.changes.values())
        #     all_mutations = [tuple(row) for row in all_mutations]
        #
        #     mutations = Counter(all_mutations)
        #     for mut in mutations.keys():
        #         if mutations[mut] >= self.n_seq/2.0:
        #             seq.sequence[int(mut[0])] = int(mut[1])
        #     return seq
        #
        # def get_seq(self, sequence_id):
        #     if sequence_id in self.changed:
        #         return self.changes[sequence_id]
        #     else:
        #         return None
        #
        # def Hamming_distance(self,simulation_settings,sample,action='mean'):
        #     HDs = []
        #     for i in sample:
        #         for j in sample:
        #             if i in self.changes.keys():
        #                 changed1 = [str(k) for k in self.changes[i]]
        #             else:
        #                 changed1 = []
        #             if j in self.changes.keys():
        #                 changed2 =  [str(k) for k in self.changes[j]]
        #             else: changed2 = []
        #             HDs.append(len(set(list(changed1)) ^ set(list(changed2))))
        #     if action == 'mean':
        #         return np.mean(HDs)
        #     elif action == 'Poisson_fit':
        #         poiss = np.mean(HDs)/(2*simulation_settings['mut_rate']*simulation_settings['seq_len'])
        #         exp = scats.poisson.pmf(range(max(HDs)+1),np.mean(HDs))*len(HDs)
        #         obs = np.histogram(HDs, bins=range(0,max(HDs)+2))[0]
        #         pval =scats.chisquare(obs,exp,ddof=len(exp)-1-len(sample)).pvalue
        #         if np.isnan(pval) or pval>0.05:
        #             return poiss
        #         else:
        #             return np.nan

if __name__ == '__main__':
    sim = Simulation(sequence=Seq(seq='ATA'))
    test = Population(sim,10)
    test.add_sequence()
    test.add_sequence([0.1,2.1])
    test.add_sequence([2.1])
    test.add_sequence([0.3,2.3])
    print test.stats()
