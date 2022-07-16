from .impactor import *

class Degrader(Impactor):
    """
    A class that represents an enzyme that degrades a substrate (or substrates).
    ...
    
    Attributes
    ----------
    enzyme : Regulator | GeneProduct
        The enzyme that degrades substrate
    substrate : Regulator | Reporter | GeneProduct | List [Regulator | Reporter | GeneProduct]
        The substrate that is degraded by enzyme
    Km : int | float | List [int | float]
        Michaelis constant
    k2 : int | float | List [int | float]
        Degradation constant to reduce enzyme-substrate complex into enzyme only
    name : str, optional
        Name of the operator displayed on the network representation
    color: str, optional
        Color displayed on the network representation

    Methods
    ----------
    degradation rate() - is calculated based on Michaelis-Menten law. Enzyme-substrate 
    complex is assumed to be in quasi-steady-state (its concenration barely changes 
    over time).
    
    """

    def __init__(self, enzyme, substrate, Km, k2, name=None, uri=None, sbol_comp=None, color='skyblue'):
        super().__init__(enzyme, name, color)
        self.enzyme = enzyme

        if type(Km)==list:
            self.km = Km
        else:
            self.km = [Km]

        if type(k2)==list:
            self.k2 = k2
        else:
            self.k2 = [k2]

        if type(substrate)==list:
            self.substrate = substrate
        else:
            self.substrate = [substrate]

    def __str__(self):
        if self.name == None:
            return 'DEGRADER'
        else: return self.name
        
    def degradation_rate(self):
        degradation_rate = []
        enzyme = []
        
        for i, substrate in enumerate(self.substrate):
            # test
            # print(f'{(self.enzyme.concentration / len(self.substrate))} split enzyme conc')

            # if enzyme has multiple substrates, enzyme is split equally between all substrates
            enzyme.append(self.enzyme.concentration / len(self.substrate))
            # calculate Vmax
            vmax = self.k2[i] * enzyme[i]
            # calculate substrate degradation rate
            substrate_change_rate = (vmax * substrate.concentration) / (self.km[i] + substrate.concentration)
            # test
            # print(f'degr_rate is {substrate_change_rate}')
            degradation_rate.append(substrate_change_rate)

        return degradation_rate