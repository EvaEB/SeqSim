"""A setuptools based setup module.
See:
https://packaging.python.org/guides/distributing-packages-using-setuptools/
https://github.com/pypa/sampleproject
"""

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

    packages=find_packages(exclude=['contrib', 'docs', 'tests']),

    python_requires='>=2.7',
    install_requires=['matplotlib==2.2','appJar','pyyaml','numpy==1.15','tqdm','scipy==1.1'],

    package_data={
        'Scenarios': ['settings_files/*'],
        'SeqSimEvo': ['simulation_settings/*']
    },

    entry_points={
        'console_scripts': [
            'SeqSimEvo=GUI:gui_run',
        ],
    },

    project_urls={  # Optional
        'Bug Reports': '',
        'Funding': '',
        'Source': '',
    },
)
