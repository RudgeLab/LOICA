import os
#from distutils.core import setup
from setuptools import setup
import subprocess

setup(name='loica',
    install_requires=[        
        'numpy', 
	'networkx',
	'pandas'
        ],
    setup_requires=[
        'numpy',
	'networkx',
	'pandas' 
        ],
    packages=['loica'],
    python_requires='>=3',
    version='v0.0'
)
