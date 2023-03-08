from .operator import *

class Sum(Operator):
    color = 'skyblue'
    shape = 's'
    def __init__(self, input, output, alpha, K, n, uri=None, sbol_comp=None):
        self.alpha = alpha
        self.K = K
        self.n = n
        self.input = input
        self.output = output
        self.uri = uri
        self.sbol_comp = sbol_comp

    def __str__(self):
        return 'SUM'
        
    def expression_rate(self, t, dt):
        expression_rate = 0
        if type(self.input)!=list:
            self.input = [self.input]
        for i,input in enumerate(self.input):
            input_regulator = input.concentration
            r = (input_regulator/self.K[i])**self.n[i]
            expression_rate += ( self.alpha[i][0] + self.alpha[i][1]*r ) / (1 + r)
        return expression_rate
