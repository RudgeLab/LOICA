class Hasty:
    color = 'orange'
    shape = 's'
    def __init__(self, input, output, alpha, sigma, tau, uri=None, sbol_comp=None):
        self.alpha = alpha
        self.sigma = sigma
        self.tau = tau
        self.input = input
        self.output = output
        self.uri = uri
        self.sbol_comp = sbol_comp

    def __str__(self):
        return 'HASTY'

    def expression_rate(self, t, dt):
        r1 = self.input[0].concentration
        r2 = self.input[1].concentration
        num = 1 + r1**2 + self.alpha * self.sigma * r1**4
        denom = (1 + r1**2 + self.sigma * r1**4) * (1 + r2**4)
        return num / denom / self.tau