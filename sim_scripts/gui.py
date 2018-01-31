from appJar import gui
import yaml
import scenarios


results = 'simulation not run yet'
final_settings = {'Organism': '','Scenario': ''}

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
    global results
    if btn == 'run':
        if final_settings['Scenario'] == 'control':
            results = scenarios.control(final_settings,plot=True)
            settings.setMessage('fasta', results)
    if btn == 'save Fasta':
        dest = settings.saveBox('save fasta as')
        with open(dest,'w') as f:
            f.write(results)

def next(btn):
    global final_settings
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
                        if type(pos_settings[i][j]) is float:
                            settings.setEntry(j,'{:f}'.format(pos_settings[i][j])) #fix for problem parsing scientific notation
                        else:
                            settings.setEntry(j,pos_settings[i][j])
                else:
                    if type(pos_settings[i]) is float:
                        settings.setEntry(i,'{:f}'.format(pos_settings[i])) #fix for problem parsing scientific notation
                    else:
                        settings.setEntry(i,pos_settings[i])

    if final_settings['Scenario'] != cur_settings['Scenario']: #Scenario has changed
        settings.setLabel('l2', 'Organism settings: {}'.format(cur_settings['Scenario']))
        with open('/home/eva/code/SeqSim/sim_scripts/settings_files/{}'.format(cur_settings['Scenario'])) as f:
            pos_settings = yaml.safe_load(f)
        for i in pos_settings:
            if type(pos_settings[i]) is dict:
                for j in pos_settings[i]:
                    if type(pos_settings[i][j]) is float:
                        settings.setEntry(j,'{:f}'.format(pos_settings[i][j])) #fix for problem parsing scientific notation
                    else:
                        settings.setEntry(j,pos_settings[i][j])
            else:
                if type(pos_settings[i]) is float:
                    settings.setEntry(i,'{:f}'.format(pos_settings[i])) #fix for problem parsing scientific notation
                else:
                    settings.setEntry(i,pos_settings[i])


    #if settings.getPagedWindowPageNumber('settings') == 3: #Scenario settings


    if settings.getPagedWindowPageNumber('settings') == 4: # Overview page
        message = ''
        for i in final_settings:
            if type( final_settings[i]) is dict:
                message+=i+'\n'
                for j in  final_settings[i]:
                    message+='\t'+ str(final_settings[i][j])+'\n'
            else:
                message+=i+'\t'+ str(final_settings[i])+'\n'
        settings.setMessage('overview_settings',message)

    for i in cur_settings:
        if i in final_settings.keys():
            final_settings[i] = yaml.load(cur_settings[i])
        elif i in final_settings['parameters'].keys():
            print cur_settings[i]
            final_settings['parameters'][i] = yaml.load(cur_settings[i])
        else:
            print i






with gui('simulation') as settings:
    settings.setResizable(canResize=True)
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

        with settings.page():
            settings.addLabel('Overview')
            with settings.scrollPane('scroll'):
                settings.addEmptyMessage('overview_settings')
            settings.addButton('run',press)

        with settings.page():
            settings.addLabel('results')
            with settings.scrollPane('scroll2'):
                settings.addMessage('fasta', results)
            settings.addButton('save Fasta', press)
