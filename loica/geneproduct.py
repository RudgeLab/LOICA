import numpy as np

class GeneProduct:
    shape = '^'
    def __init__(self, name, init_concentration=0, degradation_rate=0, uri=None, sbol_comp=None):
        self.init_concentration = init_concentration
        self.concentration = init_concentration
        self.degradation_rate = degradation_rate
        self.name = name
        self.expression_rate = 0
        self.uri = uri
        self.sbol_comp = sbol_comp


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
    color = 'lightgreen'
    def __init__(self, name, init_concentration=0, degradation_rate=0, sbol_comp=None):
        super().__init__(name, init_concentration, degradation_rate)
        self.sbol_comp = sbol_comp

class Reporter(GeneProduct):
    def __init__(self, name, init_concentration=0, degradation_rate=0, signal_id=None, color='w', sbol_comp=None):
        super().__init__(name, init_concentration, degradation_rate, sbol_comp=None)
        self.signal_id = signal_id
        self.color = color
        self.sbol_comp = sbol_comp

