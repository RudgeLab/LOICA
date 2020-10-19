class Or:
    '''
    2-input genetic OR gate
    '''
    def __init__(self, input, output, a, b, K, n, alpha=[0,1,1,0]):
        self.alpha = alpha
        self.a = a
        self.b = b
        self.K = K
        self.n = n
        self.input = input
        self.output = output
        
    def expression_rate(self, t, dt):
        input_inducer1 = self.input[0].concentration
        input_inducer2 = self.input[1].concentration
        i1 = (input_inducer1/self.K[0])**self.n[0]
        i2 = (input_inducer2/self.K[1])**self.n[1]
        expression_rate1 = ( self.a[0] + self.b[0]*i1) / (1 + i1)
        expression_rate2 = (self.a[1] + self.b[1]*i2) / (1 + i2)
        a0 = self.alpha[0]
        a1 = self.alpha[1] * expression_rate1
        a2 = self.alpha[2] * expression_rate2
        a3 = self.alpha[3] * expression_rate1 * expression_rate2
        expression_rate = a0 + a1 + a2 + a3
        return expression_rate



