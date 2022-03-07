from .operator import *

class Hasty(Operator):
    def __init__(self, input, output, alpha, sigma, tau, unit, name=None, uri=None, sbol_comp=None, color='orange'):
        super().__init__(output, name, uri, sbol_comp, color, unit)
        self.alpha = alpha
        self.sigma = sigma
        self.tau = tau
        self.input = input

    def __str__(self):
        if self.name == None:
            return 'Hasty'
        else: return self.name

    def expression_rate(self, t, dt):
        r1 = self.input[0].concentration
        r2 = self.input[1].concentration
        num = 1 + r1**2 + self.alpha * self.sigma * r1**4
        denom = (1 + r1**2 + self.sigma * r1**4) * (1 + r2**4)
        return num / denom / self.tau