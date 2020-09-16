import numpy as np

class TU:
    def output(self, inputs, dt):
        return 0

class RepressedRepressor(TU):
    def __init__(self, input, a, b, K, n):
        self.a = a
        self.b = b
        self.K = K
        self.n = n
        self.input = input
        
    def output(self, inputs, dt):
        input_repressor = inputs[self.input]
        r = (input_repressor/self.K)**self.n
        expression_rate = ( self.a + self.b*r ) / (1 + r)
        return expression_rate
