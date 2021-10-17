class Nor:
    color = 'orange'
    shape = 's'
    def __init__(self, input, output, alpha, a, b, K, n):
        self.alpha = alpha
        self.a = a
        self.b = b
        self.K = K
        self.n = n
        self.input = input
        self.output = output

    def __str__(self):
        return 'NOR'

    def expression_rate(self, t, dt):
        input_repressor1 = self.input[0].concentration
        input_repressor2 = self.input[1].concentration
        r1 = (input_repressor1/self.K[0])**self.n[0]
        r2 = (input_repressor2/self.K[1])**self.n[1]
        expression_rate1 = ( self.a[0] + self.b[0]*r1) / (1 + r1)
        expression_rate2 = (self.a[1] + self.b[1]*r2) / (1 + r2)
        a0 = self.alpha[0]
        a1 = self.alpha[1] * expression_rate1
        a2 = self.alpha[2] * expression_rate2
        a3 = self.alpha[3] * expression_rate1 * expression_rate2
        expression_rate = a0 + a1 + a2 + a3
        return expression_rate


