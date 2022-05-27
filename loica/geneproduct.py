import numpy as np

class GeneProduct:
    """
    A class that represents a gene product, protein or RNA.

    ...
    
    Attributes
    ----------
    name : str
        Name of the gene product
    init_concentration : int | float
        Initial concentration of the gene product in Molar
    degradation_rate : int | float
        Degradation rate of the gene product
    type_ : str, optional
        Molecular type of the gene product, could be 'PRO' or 'RNA'
    uri : str, optional
        SynBioHub URI
    sbol_comp : SBOL Component, optional
        SBOL Component
    """
    shape = '^'
    def __init__(self, name, init_concentration=0, degradation_rate=0, diffusion_rate=0, uri=None, sbol_comp=None, type_='PRO', color='silver'):
        self.ext_conc = 0
        self.init_concentration = init_concentration
        self.concentration = init_concentration
        self.degradation_rate = degradation_rate
        self.diffusion_rate=diffusion_rate
        self.name = name
        self.expression_rate = 0
        self.uri = uri
        self.sbol_comp = sbol_comp
        self.type_ = type_ 
        self.color = color


    def initialize(self):
        self.concentration = self.init_concentration

    def express(self, rate):
        self.expression_rate += rate


    """ 
    I need this method to account for diffusion.

    To avoid dragging Sample in here, allow input of extracellular concentration and
    determine the difference in extracellular concentration that one cell makes
    (dext_conc_dt) (might change name)
    """

    # added cells - number of cells that create this gene product. I want to conncect 
    # this to OD
    def step(self, growth_rate, dt, biomass):
        # this is how much diffused in/out of the cell
        dext_conc_dt = self.diffusion_rate*(self.concentration-self.ext_conc)
        # change of concentration within cell
        dconcdt = self.expression_rate - (self.degradation_rate + growth_rate) * self.concentration - dext_conc_dt
        self.next_concentration = self.concentration + dconcdt * dt
        self.concentration = self.next_concentration
        # then external concentration gets updated based on number of cells that produce
        # or intake this molecule
        self.next_ext_conc = self.ext_conc + dext_conc_dt * dt * biomass
        self.ext_conc = self.next_ext_conc
        self.expression_rate = 0

    def __str__(self):
        return self.name

class Regulator(GeneProduct):
    """
    Representation of a regulatory gene product.
    Child of GeneProduct.
    """
    def __init__(self, name, init_concentration=0, degradation_rate=0, diffusion_rate=0, sbol_comp=None, color='lightgreen'):
        super().__init__(name, init_concentration, degradation_rate, diffusion_rate, sbol_comp,color='lightgreen')
        self.sbol_comp = sbol_comp

class Reporter(GeneProduct):
    """
    Representation of a regulatory gene product.

    Parameters
    ----------
    signal_id : str, optional
        Flapjack ID of the signal that the reporter is associated with.
    color : str, optional
        Color of the reporter
    """
    def __init__(self, name, init_concentration=0, degradation_rate=0, diffusion_rate=0, signal_id=None, color='w', sbol_comp=None):
        super().__init__(name, init_concentration, degradation_rate, diffusion_rate, sbol_comp)
        self.signal_id = signal_id
        self.color = color
        self.sbol_comp = sbol_comp

#Producer
