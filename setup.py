# Always prefer setuptools over distutils
from setuptools import setup, find_packages
from os import path
# io.open is needed for projects that support Python 2.7
# It ensures open() defaults to text mode with universal newlines,
# and accepts an argument to specify the text encoding
# Python 3 only projects can skip this import
from io import open

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'readme.md'), encoding='utf-8') as f:
    long_description = f.read()


setup(
    name='SeqSimEvo',
    version='0.1',
    description='Simulation of sequence evolution with selection',
    long_description=long_description,
    long_description_content_type='text/markdown',

    url='',  # TODO: add

    author='Eva Bons',
    author_email='eva.bons@gmail.com',

    packages=find_packages('.',exclude=['contrib', 'docs', 'Tests']),

    python_requires='>=2.7',
    install_requires=['matplotlib','appJar','pyyaml','numpy','tqdm','scipy'],

    package_data={
        'SeqSimEvo': ['simulation_settings/*','Scenarios/settings_files/*']
    },

    entry_points={
        'console_scripts': [
            'SeqSimEvo_gui=SeqSimEvo.GUI:gui_run',
            'SeqSimEvo_multipleCompartments=SeqSimEvo.Scenarios:multiple_compartments.main',
            'SeqSimEvo_recreateDataset=SeqSimEvo.Scenarios:recreate_dataset.main',
            'SeqSimEvo_simpleSim=SeqSimEvo.Scenarios:simple_sim.main'
        ],
    },

    project_urls={
        'Bug Reports': '',
        'Funding': '',
        'Source': '',
    },
)
