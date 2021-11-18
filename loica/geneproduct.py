import numpy as np

class GeneProduct:
    """Representation of a gene product."""
    shape = '^'
    def __init__(self, name, init_concentration=0, degradation_rate=0, uri=None, sbol_comp=None, type_='PRO'):
        """Initialize a gene product.

        :param name: Name of the gene product.
        :param init_concentration: Initial concentration of the gene product.
        :param degradation_rate: Degradation rate of the gene product.
        :param uri:  SynBioHub URI of the gene product.
        :param sbol_comp: SBOL Component vinculated to the gene product.
        :param type_: Molecular type of the gene product, cpuld be 'PRO' or 'RNA'.
        """
        self.init_concentration = init_concentration
        self.concentration = init_concentration
        self.degradation_rate = degradation_rate
        self.name = name
        self.expression_rate = 0
        self.uri = uri
        self.sbol_comp = sbol_comp
        self.type_ = type_ 


    def initialize(self):
        self.concentration = self.init_concentration

    def express(self, rate):
        self.expression_rate += rate

    def step(self, growth_rate, dt):
        dconcdt = self.expression_rate - (self.degradation_rate + growth_rate) * self.concentration
        self.next_concentration = self.concentration + dconcdt * dt
        self.concentration = self.next_concentration
        self.expression_rate = 0

    def __str__(self):
        return self.name

class Regulator(GeneProduct):
    """Representation of a regulatory gene product."""
    color = 'lightgreen'
    def __init__(self, name, init_concentration=0, degradation_rate=0, sbol_comp=None):
        super().__init__(name, init_concentration, degradation_rate)
        self.sbol_comp = sbol_comp

class Reporter(GeneProduct):
    """Representation of a regulatory gene product."""
    def __init__(self, name, init_concentration=0, degradation_rate=0, signal_id=None, color='w', sbol_comp=None):
        """Initialize Reporter
        
        :param signal_id: Flapjack ID of the signal that the reporter is associated with.
        :param color: Color of the reporter.
        """
        super().__init__(name, init_concentration, degradation_rate, sbol_comp=None)
        self.signal_id = signal_id
        self.color = color
        self.sbol_comp = sbol_comp

#Producer
