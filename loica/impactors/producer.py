from .impactor import *

class Producer(Impactor):
    """
    A class that represents an enzyme that synthesizes a product.

    Assuming that the concentration of substrate used to synthesize product is high 
    (increasing to infinity), we can say that all enzymes are engaged in 
    enzyme-substrate complex. Therefore, only k2 constant is used here.
    ...
    
    Attributes
    ----------
    enzyme : Regulator | GeneProduct
        The enzyme that degrades substrate
    product : Regulator | Reporter | GeneProduct
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
            return 'DEGRADER'
        else: return self.name
        
    def production_rate(self):
        return self.enzyme.concentration * self.k2
        