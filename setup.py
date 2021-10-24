import os
#from distutils.core import setup
from setuptools import setup
import subprocess

setup(name='loica',
    install_requires=[        
        'numpy', 
	'networkx',
	'pandas',
    'sbol3',
    'tyto'
        ],
    setup_requires=[
        'numpy',
	'networkx',
	'pandas',
    'sbol3',
    'tyto' 
        ],
    packages=['loica'],
    python_requires='>=3',
    version='v0.0'
)
