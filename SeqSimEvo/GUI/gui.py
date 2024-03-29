import os
from collections import OrderedDict
import pkg_resources

import matplotlib
matplotlib.use('TkAgg')
from appJar import gui
import yaml

from . import scenarios
import fasta_tools


def select():
    selection = app.getAllOptionBoxes()

    #display default settings organism
    global sim_settings_path
    sim_settings_path = pkg_resources.resource_filename('SeqSimEvo','simulation_settings/')
    global scenario_settings_path
    scenario_settings_path = pkg_resources.resource_filename('SeqSimEvo','Scenarios/settings_files/')
    with open(sim_settings_path+selection['Organism']) as f:
        OrgSettings = f.read()
    app.clearTextArea('OrgSet')
    app.setTextArea('OrgSet', OrgSettings)

    #display default settings Scenario
    with open(scenario_settings_path+selection['Scenario']) as f:
        ScenSettings = f.read()
    app.clearTextArea('ScenSet')
    app.setTextArea('ScenSet', ScenSettings)

def update_fasta(fasta):
    app.clearTextArea('fasta_text')
    app.setTextArea('fasta_text',fasta)

def settings_press(btn):
    settings = app.getAllTextAreas()
    major_settings = app.getAllOptionBoxes()
    if btn == 'save organism settings':
        #specify file location
        filename = app.saveBox(title='Save organism settings',
                               dirName=sim_settings_path,
                               fileExt="")
        #save file
        with open(filename, 'w') as f:
            f.write(settings['OrgSet'])
    elif btn == 'save scenario settings':
        #specify file location
        filename = app.saveBox(title='Save scenario settings',
                               dirName=scenario_settings_path,
                               fileExt="")
        #save file
        with open(filename, 'w') as f:
            f.write(settings['ScenSet'])
    elif btn == 'run':
        fasta = app.threadCallback(scenarios.run,update_fasta,major_settings['Scenario'],
                                  yaml.safe_load(settings['ScenSet']),
                                  yaml.safe_load(settings['OrgSet']))
    elif btn == 'save fasta':
        #specify file location
        filename = app.saveBox(title='Save fasta',fileExt=".fasta")
        #save file
        with open(filename, 'w') as f:
            f.write(settings['fasta_text'])

    elif btn == 'view highlighter':
        app.thread(fasta_tools.highlighter(fasta_string=settings['fasta_text'],consensus=0))

    elif btn == 'view nj tree':
        app.thread(fasta_tools.njTree(fasta_string=settings['fasta_text']))

def gui_run():
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    global app
    with gui('SeqSim') as app:
        dir_path = os.path.dirname(os.path.realpath(__file__))
        dir_path_up = os.sep.join(dir_path.split(os.sep)[:-1])
        sim_settings_path = pkg_resources.resource_filename('SeqSimEvo','simulation_settings/')
        scenario_settings_path = pkg_resources.resource_filename('SeqSimEvo','Scenarios/settings_files/')
        Organism_options = os.listdir(sim_settings_path)
        app.addLabelOptionBox('Organism',Organism_options)
        Scenario_options = os.listdir(scenario_settings_path)
        app.addLabelOptionBox('Scenario',Scenario_options)

        app.addButton('select',select)

        app.addLabel('Org','Organism settings')
        app.setLabelBg('Org','grey')
        app.addTextArea('OrgSet')
        app.getTextAreaWidget('OrgSet').config(font="Courier 20")
        app.addButton('save organism settings',settings_press)

        app.addLabel('Scen','Scenario settings')
        app.setLabelBg('Scen','grey')
        app.addTextArea('ScenSet')
        app.getTextAreaWidget('ScenSet').config(font="Courier 20")
        app.addButtons(['save scenario settings','run'],settings_press)

        app.addLabel('fasta','simulation output')
        app.setLabelBg('fasta','grey')
        app.addTextArea('fasta_text')
        app.getTextAreaWidget('fasta_text').config(font="Courier 20")
        app.addButtons(['save fasta','view highlighter','view nj tree'],settings_press)

if __name__ == '__main__':
    gui_run()
