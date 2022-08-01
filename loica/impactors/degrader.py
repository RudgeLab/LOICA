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
    degradation rate() - is calculated based on Michaelis-Menten law and equation 4 
    (competing substrates kinetics) from https://doi.org/10.1016/j.febslet.2013.06.025. 
    Enzyme-substrate complex is assumed to be in quasi-steady-state (its concenration barely changes 
    over time).
    
    """

    def __init__(self, enzyme, substrate, Km, k2, mode=None, name=None, uri=None, sbol_comp=None, color='skyblue'):
        super().__init__(enzyme, name, color)
        self.enzyme = enzyme

        if type(substrate)==list:
            self.substrate = substrate
        else:
            self.substrate = [substrate]

        if type(Km)==list:
            self.km = Km
            if len(self.km)>len(self.substrate) or (len(self.km)<len(self.substrate) and len(self.km)!=1):
                print('There should be either the same number of Km constants as substrates, or one Km value for all')
        else:
            self.km = []
            for s in self.substrate:
                self.km.append(Km)

        if type(k2)==list:
            self.k2 = k2
            if len(self.km)>len(self.substrate) or (len(self.k2)<len(self.substrate) and len(self.k2)!=1):
                print('There should be either the same number of k2 parameters as substrates, or one k2 value for all')
        else:
            self.k2 = []
            for s in self.substrate:
                self.k2.append(k2)

        self.mode=mode


    def __str__(self):
        if self.name == None:
            return 'DEGRADER'
        else: return self.name
        
    def degradation_rate(self):
        degradation_rate = []

        for i, substrate in enumerate(self.substrate):
            x = 1
            # calculate Vmax
            vmax = self.k2[i] * self.enzyme.concentration

            # calculate substrate degradation rate
            for ii, s in enumerate(self.substrate):
                if i != ii:
                    if self.mode=='papers':
                        substrate_change_rate = (vmax * substrate.concentration) / (self.km[i] + substrate.concentration + s.concentration)
                    elif not self.mode:
                        x+= s.concentration / self.km[ii]
    
            if not self.mode:
                substrate_change_rate = (vmax * substrate.concentration) / (self.km[i]*x + substrate.concentration)
            # test
            # print(f'degr_rate is {substrate_change_rate}')
            degradation_rate.append(substrate_change_rate)

        return degradation_rate