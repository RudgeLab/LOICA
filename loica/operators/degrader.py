from .operator import *
import numpy as np
from scipy.optimize import least_squares
from .receiver import *

class Degrader(Operator):
    """
    A class that represents an enzyme that degrades a substrate (or substrates).
    ...
    
    Attributes
    ----------
    input : Regulator | GeneProduct
        The enzyme that degrades substrate
    output : Regulator | Reporter | GeneProduct | List [Regulator | Reporter | GeneProduct]
        The substrate that is degraded by enzyme
    k1 : int | float | List [int | float]
        Binding constant to produce enzyme-substrate complex
    k1r : int | float | List [int | float]
        Dissociation constant to reduce the enzyme-substrate complex
    k2 : int | float | List [int | float]
        Degradation constant to reduce enzyme-substrate complex into enzyme only
    name : str, optional
        Name of the operator displayed on the network representation
    color: str, optional
        Color displayed on the network representation

    """

    def __init__(self, input, output, k1, k1r, k2, name=None, uri=None, sbol_comp=None, color='skyblue'):
        super().__init__(output, name, uri, sbol_comp, color)
        self.input = input

        if type(k1)==list:
            self.k1 = k1
        else:
            self.k1 = [k1]

        if type(k1r)==list:
            self.k1r = k1r
        else:
            self.k1r = [k1r]

        if type(k2)==list:
            self.k2 = k2
        else:
            self.k2 = [k2]

        if type(output)==list:
            self.output = output
        else:
            self.output = [output]

        self.es_complex = []
        for o in self.output:
            self.es_complex.append(0)

    def __str__(self):
        if self.name == None:
            return 'DEGRADER'
        else: return self.name
        
    def degradation_rate(self, dt):
        degradation_rate = []
        enzyme = []
        
        for i, output in enumerate(self.output):
            # test
            print(f'{(self.input.concentration / len(self.output) - self.es_complex[i])} split enzyme conc')
            # if enzyme has multiple substrates, enzyme is split equally between all substrates (account for enzymes in complex)
            enzyme.append(self.input.concentration / len(self.output) - self.es_complex[i])

            # calculate substrate degradation rate
            substrate_change_rate = -self.k1 * enzyme[i] * output.concentration + self.k1r * self.es_complex[i] - self.k2 * self.es_complex[i]
            # test
            print(f'degr_rate is {substrate_change_rate}')
            degradation_rate.append(-substrate_change_rate)

            # update enzyme-substrate complex concentrations
            es_change = self.k1 * enzyme[i] * output.concentration - self.k1r * self.es_complex[i] - self.k2 * self.es_complex[i]
            # enzyme_change = -self.k1 * enzyme[i] * output.concentration + self.k1r * self.es_complex[i] + self.k2 * self.es_complex[i]
            new_es = self.es_complex[i] + es_change * dt
            self.es_complex[i] = new_es

        return degradation_rate