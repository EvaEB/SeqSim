 import matplotlib
matplotlib.use('TkAgg')
from appJar import gui
import yaml
import os
import scenarios
from collections import OrderedDict
import fasta_tools
from IPython.utils.capture import capture_output


def select():
    selection = app.getAllOptionBoxes()

    #display default settings organism
    with open('../seq_sim/simulation_settings/'+selection['Organism']) as f:
        OrgSettings = f.read()
    app.clearTextArea('OrgSet')
    app.setTextArea('OrgSet', OrgSettings)

    #display default settings Scenario
    with open('settings_files/'+selection['Scenario']) as f:
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
                               dirName='../seq_sim/simulation_settings',
                               fileExt="")
        #save file
        with open(filename, 'w') as f:
            f.write(settings['OrgSet'])
    elif btn == 'save scenario settings':
        #specify file location
        filename = app.saveBox(title='Save scenario settings',
                               dirName='settings_files/',
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

os.chdir(os.path.dirname(os.path.abspath(__file__)))
with gui('SeqSim') as app:
    Organism_options = os.listdir('../seq_sim/simulation_settings/')
    app.addLabelOptionBox('Organism',Organism_options)

    Scenario_options = os.listdir('settings_files/')
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
