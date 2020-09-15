class TU:
    def output(self, inputs, dt):
        return 0

class RepressedRepressor(TU):
    def __init__(self, a, b, K, n):
        self.a = a
        self.b = b
        self.K = K
        self.n = n
        
    def output(self, input_repressor, dt):
        r = (input_repressor/self.K)**self.n
        expression_rate = ( self.a + self.b*r ) / (1 + r)
        return expression_rate
