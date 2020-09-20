
class Receiver:
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

