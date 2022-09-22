from .impactor import *

class Producer(Impactor):
    """
    A class that represents an enzyme active site that synthesizes a product.

    Assuming that the concentration of substrate used to synthesize product is high 
    (increasing to infinity), Km+[S]=[S]. This leads to Km+[S] and [S] cancelling each 
    other in Michaelis-Menten equation, so only k2 constant is needed to calculate 
    synthesis rate.
    ...
    
    Attributes
    ----------
    enzyme : Regulator 
        The enzyme that degrades substrate
    product : Regulator | Reporter 
        The product synthesized by enzyme
    k2 : int | float
       Synthesis constant to synthesize product from enzyme-substrate complex
    name : str, optional
        Name of the operator displayed on the network representation
    color: str, optional
        Color displayed on the network representation

    """

    def __init__(self, enzyme, product, k2, name=None, uri=None, sbol_comp=None, color='skyblue'):
        super().__init__(enzyme, name, color)
        self.enzyme = enzyme
        self.k2 = k2
        self.product = product

    def __str__(self):
        if self.name == None:
            return 'PRODUCER'
        else: return self.name
        
    def production_rate(self):
        return self.enzyme.concentration * self.k2
        