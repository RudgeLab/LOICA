from setuptools import setup, find_packages


with open('README.md', 'r') as ld:
    long_description = ld.read()

setup(name='loica',
    version='1.0.5',
    author='Gonzalo Vidal',
    author_email='gsvidal@uc.cl',
    description='LOICA: Logical Operators for Integrated Cell Algorithms',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/RudgeLab/LOICA',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        ],
    install_requires=[        
        'numpy', 
	    'scipy',
	    'networkx',
	    'pandas',
    	'sbol3',
    	'tyto',
   	    'tdqm'
        ],
    setup_requires=[
        'numpy',
	    'scipy',
	    'networkx',
	    'pandas',
    	'sbol3',
    	'tyto',
    	'tdqm'
        ],
    packages=find_packages(),
    python_requires='>=3'
)
