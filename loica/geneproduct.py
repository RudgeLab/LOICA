import numpy as np

class GeneProduct:
    def __init__(self, name, init_concentration=0, degradation_rate=0):
        self.concentration = init_concentration
        self.degradation_rate = degradation_rate
        self.name = name
        self.expression_rate = 0

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
    def __init__(self, name, init_concentration=0, degradation_rate=0):
        super().__init__(name, init_concentration, degradation_rate)

class Reporter(GeneProduct):
    def __init__(self, name, init_concentration=0, degradation_rate=0, signal_id=None):
        super().__init__(name, init_concentration, degradation_rate)
        self.signal_id = signal_id


