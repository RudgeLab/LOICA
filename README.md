# LOICA
Logical Operators for Intelligent Cell Algorithms

## High level idea
Program cells with Python

Easy progamation of genetic network models

Generation of synthetic data

Generation of Human/Machine readable protocols for DNA assembly

Communication with FlapJack

Use and output SBOL/SBML

Use all sorts of cellular computation


## Usage

```python

r1 = Repressor(input=0, output=1, a=1e4, e=0, n=2, dt=0.05) 
r2 = Repressor(input=0, output=1, a=1e4, e=0, n=2, dt=0.05) 
r3 = Repressor(input=0, output=1, a=1e4, e=0, n=2, dt=0.05) 
repressilator = Repressilator(proteins= [5,0,0], r1=r1, r2=r2, r3=r3)
metabolism = Metabolism(gamma = 0.1, mu=1, dt=0.05)
cell = Cell(metabolism=metabolism, genetic_network=repressilator)

```
Ploting results

```python
cell.run(2000)
memory_np = np.array(cell.memory)
ap1 = memory_np[:,0]
ap2 = memory_np[:,1]
ap3 = memory_np[:,2]

t = np.arange(2000) * 0.05
plt.plot(t, ap1, 'r')
plt.plot(t, ap2, 'g')
plt.plot(t, ap3, 'b')

plt.xlabel('Time')
```
Output
<img src="https://github.com/SynBioUC/LOICA/blob/master/images/time_dynamics.png" height="300" />


Nice!
