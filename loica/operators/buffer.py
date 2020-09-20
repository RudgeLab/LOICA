
class Buffer:
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


