from __future__ import print_function
import SeqSimEvo
#from SeqSimEvo import binary_fission as SeqSimEvo

import pytest
import numpy as np
from copy import deepcopy, copy

#tests written for SeqSimEvo.Seq.__init__
#                               .__str__
#                               .__getitem__
#                               .translate_sequence
#                               .generate_seq


#__init__
def test_Seq_init_sequence_correct_length():
    sequence = SeqSimEvo.Seq()

    assert(len(sequence.sequence) == sequence.len)

def test_Seq_init_creation():
    '''test for presence of the necessary variables'''
    sequence = SeqSimEvo.Seq()

    sequence.translation
    sequence.len
    sequence.sequence

#__str__
def test_Seq_str():
    sequence = 'ATGCTGC'
    assert(str(SeqSimEvo.Seq(seq=sequence))==sequence)

#__getitem__
def test_Seq_getitem_CORRECT():
    orignal_sequence = 'ATGC'
    sequence = SeqSimEvo.Seq(seq=orignal_sequence)
    assert sequence[0] == 'A', 'returned {} when \'A\' was expected'.format(sequence[0])
    assert sequence[1] == 'T', 'returned {} when \'T\' was expected'.format(sequence[0])
    assert sequence[2] == 'G', 'returned {} when \'G\' was expected'.format(sequence[0])
    assert sequence[3] == 'C', 'returned {} when \'C\' was expected'.format(sequence[0])

def test_Seq_getitem_FAILS():
    sequence = SeqSimEvo.Seq()
    with pytest.raises(Exception) as e_info:
        sequence[len(sequence)+1]

#translate_sequence
def test_Seq_translate_sequence_CORRECT():
    orignal_sequence = 'ATCGTTCG'
    expected_translation = [0,2,3,1,2,2,3,1]
    result = SeqSimEvo.Seq().translate_sequence(orignal_sequence)
    assert list(result) == expected_translation, '{} != {}'.format(result,expected_translation)

#generate_seq
def test_Seq_generate_seq_CORRECT():
    seq_len = 123
    sim = SeqSimEvo.Seq()
    seq = sim.generate_seq(seq_len)
    assert len(sim.sequence) == seq_len, 'sequence length should be 123, is {}'.format(len(seq))


    seq = sim.generate_seq(5,base_dist=[1,1,1,1])
    assert sim.sequence == [0,0,0,0,0]

def test_Seq_generate_seq_FAILS_negLength():
    seq_len = -123
    sim = SeqSimEvo.Seq()
    with pytest.raises(Exception) as e_info:
        seq = sim.generate_seq(seq_len)

def test_Seq_generate_seq_FAILS_nonIntLenght():
    seq_len = 5.5
    sim = SeqSimEvo.Seq()
    with pytest.raises(Exception) as e_info:
        seq = sim.generate_seq(seq_len)

# tests written for SeqSimEvo.Population.__init__
#                                       .__len__
#                                       .__str__
#                                       .add_change
#                                       .add_sequence
#                                       .consensus_sequence
#                                       .copy
#                                       .delete_sequence
#                                       .get_base
#                                       .get_sample
#                                       .get_seq
#                                       .print_sample
#                                       .sample_to_string
#                                       .stats
#                                       .to_fasta
#                                       .Hamming_distance

#__init__
def test_Population_init_CORRECT():
    sim = SeqSimEvo.Simulation()
    pop = SeqSimEvo.Population(sim)
    pop.changed
    pop.changes
    pop.sim
    pop.n_seq


#__len__
def test_Population_len_CORRECT():
    sim = SeqSimEvo.Simulation(n_seq_init=100)
    pop = SeqSimEvo.Population(sim)
    assert len(pop) == 100
    assert pop.n_seq == 100

#__str__
def test_Population_str_CORRECT():
    #tests that a string gets returned, not the contents
    sim = SeqSimEvo.Simulation(n_seq_init=100)
    pop = SeqSimEvo.Population(sim)
    assert type(str(pop)) == type(''), "str() doesn't return a string, but {}".format(type(str(sim)))

#add_change
def test_Population_add_change_CORRECT():
    sim = SeqSimEvo.Simulation(n_seq_init=100,seq_len=200)
    pop = SeqSimEvo.Population(sim)
    pop.add_change(50,10,0)
    assert 50 in pop.changed, 'sequence 50 not in changed {}'.format(pop.changed)
    assert [10, 0] in pop.get_seq(50), '{}'.format(pop.changes[50])

def test_Population_add_change_FAILS_outofrange():
    sim = SeqSimEvo.Simulation(n_seq_init=100,seq_len=100)
    pop = SeqSimEvo.Population(sim)

    with pytest.raises(Exception) as e_info:
        pop.add_change(150,10,0)

    with pytest.raises(Exception) as e_info:
        pop.add_change(50,210,0)

# add_sequence
def test_Population_add_sequence_CORRECT():
    sim = SeqSimEvo.Simulation(n_seq_init=100,seq_len=100)
    pop = SeqSimEvo.Population(sim)

    new_change = [[10,1],[20,2]]
    new = pop.add_sequence(changes=new_change)
    assert new == 100, 'new sequence has wrong ID assigned. is {}, should be 100'.format(new)
    assert ([10,1] in pop.changes[new]) & ([20,2] in pop.changes[new]), 'changes weren\'t added correctly. is {}, should be {}'.format(pop.changes[new],new_change)

#consensus_sequence
def test_Population_consensus_sequence_CORRECT():
    seq = SeqSimEvo.Seq(seq='AAAA')
    sim = SeqSimEvo.Simulation(n_seq_init=1,sequence=seq)
    pop = SeqSimEvo.Population(sim)
    for i in range(10):
        pop.add_sequence(changes=[[1,2]])
    assert str(pop.consensus_sequence()) == 'ATAA', 'consensus sequence wrong. is {}, should be {}'.format(pop.consensus_sequence(),'ATAA')

#copy
def test_Population_copy_CORRECT():
    seq = SeqSimEvo.Seq(seq='AAAAAAAAAA')
    sim = SeqSimEvo.Simulation(n_seq_init=1,sequence=seq)
    pop = SeqSimEvo.Population(sim)

    for i in range(10):
        pop.add_sequence(changes=[[i,2]])

    pop_copied = pop.copy()

    for i in range(10):
        pop_copied.add_sequence(changes=[[i,3]])

    assert len(pop_copied) != len(pop)
    assert pop_copied.changes != pop.changes

#delete_sequence
def test_Population_delete_sequence_CORRECT():
    seq = SeqSimEvo.Seq(seq='AAAAAAAAAA')
    sim = SeqSimEvo.Simulation(n_seq_init=1,sequence=seq)
    pop = SeqSimEvo.Population(sim)

    id1 = pop.add_sequence(changes=[[2,2]])
    id2 = pop.add_sequence(changes=[[3,3]])

    pop.delete_sequence(id1)

    assert [3,3] in pop.get_seq(id1), '{}'.format(id2)
    assert pop.n_seq == 2

#get_base
def test_Population_get_base_CORRECT():
    seq = SeqSimEvo.Seq(seq='AAAAAAAAAA')
    sim = SeqSimEvo.Simulation(n_seq_init=1,sequence=seq)
    pop = SeqSimEvo.Population(sim)

    assert pop.get_base(0,1) == 0

#get_sample
def test_Population_get_sample_CORRECT():
    seq = SeqSimEvo.Seq(seq='AAAAAAAAAA')
    sim = SeqSimEvo.Simulation(n_seq_init=100,sequence=seq)
    pop = SeqSimEvo.Population(sim)

    sample = pop.get_sample(10)
    assert len(sample) == 10, 'sample returns wrong length'
    assert len(set(sample)) == 10,'sample returns double seqIDS'

    sample = pop.get_sample(500)
    assert len(sample) == 100, 'sample returns wrong length when oversampling'
    assert len(set(sample)) == 100,'sample returns double seqIDS when oversampling'

#get_seq
def test_Population_get_seq_CORRECT():
    seq = SeqSimEvo.Seq(seq='AAAAAAAAAA')
    sim = SeqSimEvo.Simulation(n_seq_init=100,sequence=seq)
    pop = SeqSimEvo.Population(sim)

    id = pop.add_sequence(changes=[[5,2]])

    assert [5,2] in pop.get_seq(id)
    assert pop.get_seq(0) is None


    with pytest.raises(Exception) as e_info:
        pop.get_seq(200)

#print_sample
def test_Population_print_sample_correct(capfd):
    seq = SeqSimEvo.Seq(seq='AAAAAAAAAA')
    sim = SeqSimEvo.Simulation(n_seq_init=100,sequence=seq)
    pop = SeqSimEvo.Population(sim,changes={1:[[5,2],[3,1]]},changed=[1])

    pop.print_sample([1])
    out, err = capfd.readouterr()

    assert ('A-5-T\t1' in out) and ('A-3-G\t1' in out), 'printed output:\n{}\nis incomplete. should contain A-5-T\t1 and A-3-G\t1'.format(out, '')

    with pytest.raises(Exception) as e_info:
        pop.print_sample([200])

#sample_to_string
def test_Population_sample_to_string_correct():
    seq = SeqSimEvo.Seq(seq='AAAAAAAAAA')
    sim = SeqSimEvo.Simulation(n_seq_init=100,sequence=seq)
    pop = SeqSimEvo.Population(sim,changes={1:[[5,2],[3,1]]},changed=[1])

    string = pop.sample_to_string([1])
    assert ('A-5-T\t1' in string) and ('A-3-G\t1' in string), 'returned string:\n{}\nis incomplete. should contain A-5-T\t1 and A-3-G\t1'.format(string, '')

#stats
def test_Population_stats_correct():
    seq = SeqSimEvo.Seq(seq='AAAAAAAAAA')
    sim = SeqSimEvo.Simulation(n_seq_init=100,sequence=seq)
    pop = SeqSimEvo.Population(sim,changes={1:[[5,2],[3,1]]},changed=[1])
    for i in range(60):
        pop.add_change(i,6,3)
    pop.add_sequence(changes=[[5,2],[1,1]])

    statistics = pop.stats()

    assert statistics['n_seq'] == 101, 'N_seq statistic wrong'
    assert statistics['unmutated'] == 40, 'number unmutated sequences wrong'
    assert statistics['total_mutations'] == 64, 'total number of mutations is wrong'
    assert statistics['unique_mutations'] == 4, 'number of unique mutations is wrong'
    assert statistics['majority_mutations'] == 1, 'number of majority mutations is wrong'
    assert statistics['max_fraction'] == 60.0/101, 'maximum fraction is wrong'
    assert statistics['GA_rate'] ==  0, 'fraction G-A is wrong'

#to_fasta
def test_Population_to_fasta_correct():
    seq = SeqSimEvo.Seq(seq='AAAAAAAAAA')
    sim = SeqSimEvo.Simulation(n_seq_init=3,sequence=seq)
    pop = SeqSimEvo.Population(sim,changes={1:[[5,2],[3,1]]},changed=[1])

    fasta = '>0\nAAAAAAAAAA\n>1\nAAAGATAAAA\n>2\nAAAAAAAAAA\n'
    assert pop.to_fasta([0,1,2]) == fasta, 'fasta incorrect'

#Hamming_distance
def test_Population_hamming_distance_correct():
    seq = SeqSimEvo.Seq(seq='AAAAAAAAAA')
    sim = SeqSimEvo.Simulation(n_seq_init=10,sequence=seq)
    pop = SeqSimEvo.Population(sim,changes={1:[[5,2],[3,1]],
                                            2:[[5,1],[3,1]],
                                            3:[[4,2],[3,1]],
                                            4:[[4,2],[3,1]],
                                            5:[[4,2],[3,1],[6,1]]},
                                    changed=[1,2,3,4,5])
    #action = mean
    assert abs(pop.Hamming_distance(range(10)) - 1.6) < 0.05, 'mean hamming distance wrong'
    #action = poisson_fit
    assert pop.Hamming_distance(range(10),action = 'Poisson_fit') is np.nan

    pop = SeqSimEvo.Population(sim)
    assert pop.Hamming_distance(range(10),action = 'Poisson_fit') == 0

#                           .Simulation.__init__
#                                      .__str__
#                                      .copy
#                                      .get_fitness_effect
#                                      .get_nr_offspring
#                                      .mutate_seq
#                                      .new_generation
#                                      .new_mutations_per_seq

#__init__
def test_Simulation_init_Correct():
    sim = SeqSimEvo.Simulation(model='lognormal',
                               parameters={'fl':0.3,'mu':0.6,'sigma':0.4},
                               mut_rate=0.001,
                               subs_matrix=[[0  ,0.1,0.6,1],
                                            [0.2,0.2,0.6,1],
                                            [0.4,0.5,0.5,1],
                                            [0.3,0.6,1  ,1]],
                               seq_len = 10,
                               basedist = [0,0.5,0.8,1],
                               R0 = 6.5,
                               max_pop = 10000,
                               name = 'testingSim')

    sim.settings
    assert sim.mut_rate == 0.001
    assert len(sim.sequence) == 10

    assert sim.settings['basedist'] == [0,0.5,0.8,1]
    assert 'sequence' not in sim.settings.keys()

    assert sim.sequence.sequence.count(0) == 0
    assert sim.gen == 0
    assert sim.average_fitness == 1
    assert sim.effective_pop == 1
    assert sim.n_seq == 1

    assert sim.current_gen.changes == {}
    sim.mutations_per_seq


def test_Simulation_get_fitness_effect():
    seq = SeqSimEvo.Seq(seq='AAAAAAAAAA')
    sim = SeqSimEvo.Simulation(n_seq_init=3,sequence=seq)

    assert sim.get_fitness_effect(0,0) == 1

    with pytest.raises(Exception) as e_info:
        sim.get_fitness_effect(50,0) == 1

    with pytest.raises(Exception) as e_info:
        sim.get_fitness_effect(5,10) == 1


#copy
def test_Simulation_copy():
    seq = SeqSimEvo.Seq(seq='AAAAAAAAAA')
    sim = SeqSimEvo.Simulation(n_seq_init=3,sequence=seq)

    sim_copy = sim.copy('sim_test_copy')

    sim.new_generation()

    assert sim.gen != sim_copy.gen

#get_nr_offspring
def test_Simulation_get_nr_offspring():
    seq = SeqSimEvo.Seq(seq='AAAAAAAAAA')
    sim = SeqSimEvo.Simulation(n_seq_init=3,sequence=seq)
    offspring = sim.get_nr_offspring(0)
    assert offspring >=0

    offspring,fitness = sim.get_nr_offspring(0,return_fitness=True)
    assert fitness == 1

#mutate_seq
def test_Simulation_mutate_seq():
    seq = SeqSimEvo.Seq(seq='AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA')
    sim = SeqSimEvo.Simulation(n_seq_init=3,sequence=seq)
    next_gen = SeqSimEvo.Population(simulation=sim)

    sim.mutations_per_seq = iter([3])
    sim.mutate_seq(next_gen,0,0)

    assert next_gen.stats()['total_mutations'] == 3

def test_Simulation_new_generation():
    seq = SeqSimEvo.Seq(seq_len=5000)
    sim = SeqSimEvo.Simulation(n_seq_init=10000,sequence=seq,mut_rate=0.001)
    sim.effective_pop = -1
    sim.average_fitness = -1
    sim.n_seq = -1
    old_gen = deepcopy(sim.current_gen.changes)

    sim.new_generation()

    assert sim.gen == 1
    assert sim.current_gen.changes != old_gen, 'this test can fail by chance. Be worried if it keeps failing.'
    assert sim.effective_pop != -1
    assert sim.average_fitness != -1
    assert sim.n_seq != -1
