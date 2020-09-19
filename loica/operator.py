import numpy as np

class Operator:
    def output(self, inputs, dt):
        return 0


class Source(Operator):
    def __init__(self, output, rate):
        self.rate = rate
        self.output = output

    def expression_rate(self):
        return self.rate

class Not(Operator):
    def __init__(self, input, output, a, b, K, n):
        self.a = a
        self.b = b
        self.K = K
        self.n = n
        self.input = input
        self.output = output
        
    def expression_rate(self):
        input_repressor = self.input.concentration
        r = (input_repressor/self.K)**self.n
        expression_rate = ( self.a + self.b*r ) / (1 + r)
        return expression_rate


class Buffer(Operator):
    def __init__(self, input, output, a, b, K, n):
        self.a = a
        self.b = b
        self.K = K
        self.n = n
        self.input = input
        self.output = output
        
    def expression_rate(self):
        input_inducer = self.input.concentration
        i = (input_inducer/self.K)**self.n
        expression_rate = ( self.a + self.b*i ) / (1 + i)
        return expression_rate

class Receiver(Operator):
    def __init__(self, receptor, inducer, output, a, b, K, n):
        self.a = a
        self.b = b
        self.K = K
        self.n = n
        self.receptor = 0
        self.inducer = inducer
        self.receptor = receptor
        self.output = output
        
    def expression_rate(self):
        receptor = self.receptor.concentration
        inducer = self.inducer.concentration
        complex = receptor * inducer
        i = (complex/self.K)**self.n
        expression_rate = ( self.a + self.b*i ) / (1 + i)
        return expression_rate

