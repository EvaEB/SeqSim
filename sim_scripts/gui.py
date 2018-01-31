from appJar import gui
import yaml


global final_settings

def form_from_file(fname,screen):
    with open(fname) as f:
        pos_settings = yaml.safe_load(f)

    for i in pos_settings:
        final_settings[i] = pos_settings[i]
        if type(pos_settings[i]) is dict:
            screen.addLabel('--')
            screen.addLabel(i,i)
            for j in pos_settings[i]:
                screen.addLabelEntry(j)
                screen.setEntry(j,pos_settings[i][j])
            screen.addLabel('---')
        else:
            screen.addLabelEntry(i)
            screen.setEntry(i,pos_settings[i])

def press(btn):
    with settings.page(windowTitle="settings", pageNumber=2):
        settings.addLabel("newLab", "New Label")

def next(btn):
    cur_settings = settings.getAllOptionBoxes()
    entries = settings.getAllEntries()
    for i in entries:
        cur_settings[i] = entries[i]
    if settings.getPagedWindowPageNumber('settings') == 2: #organism settings
        if final_settings['Organism'] != cur_settings['Organism']: #organism has changed
            settings.setLabel('l2', 'Organism settings: {}'.format(cur_settings['Organism']))
            with open('/home/eva/code/SeqSim/seq_sim/simulation_settings/{}'.format(cur_settings['Organism'])) as f:
                pos_settings = yaml.safe_load(f)
            for i in pos_settings:
                if type(pos_settings[i]) is dict:
                    for j in pos_settings[i]:
                        settings.setEntry(j,pos_settings[i][j])
                else:
                    settings.setEntry(i,pos_settings[i])

    for i in cur_settings:
        if i in final_settings.keys():
            final_settings[i] = cur_settings[i]
        elif i in final_settings['parameters'].keys():
            final_settings['parameters'][i] = cur_settings[i]
        else:
            print i





final_settings = {'Organism': '','Scenario': ''}
with gui('simulation') as settings:
    with settings.pagedWindow('settings'):
        settings.setPagedWindowFunction('settings', next)
        with settings.page():
            settings.addLabel('Basic selection')
            settings.addLabelOptionBox('Scenario',['control','-skyline','-migration','-root'])
            settings.addLabelOptionBox('Organism', ['HIV','phix174'])

        with settings.page():
            settings.addLabel('l2','Organism Settings')
            form_from_file('/home/eva/code/SeqSim/seq_sim/simulation_settings/template',settings)

        with settings.page():
            settings.addLabel('l4','Scenario Settings')
            form_from_file('/home/eva/code/SeqSim/sim_scripts/settings_files/template',settings)

for i in final_settings:
    print i, '\t', final_settings[i]
