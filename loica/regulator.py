import numpy as np

class Regulator:
    def __init__(self, name, init_concentration=0, degradation_rate=0):
        self.concentration = init_concentration
        self.degradation_rate = degradation_rate
        self.name = name

    def step(self, expression_rates, growth_rate, dt):
        expression_rate = np.sum(expression_rates)
        dconcdt = expression_rate - (self.degradation_rate + growth_rate) * self.concentration
        self.next_concentration = self.concentration + dconcdt * dt

    def update(self):
        self.concentration = self.next_concentration

    def __str__(self):
        return self.name
