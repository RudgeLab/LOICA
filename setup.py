import os
#from distutils.core import setup
from setuptools import setup
import subprocess

setup(name='loica',
    install_requires=[        
        'numpy', 
        ],
    setup_requires=[
        'numpy', 
        ],
    packages=['loica'],
    python_requires='>=3',
    version='v0.0'
)