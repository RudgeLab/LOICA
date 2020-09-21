class GeneticNetwork():
    def __init__(self, vector=None):
        self.operators = []
        self.regulators = []
        self.reporters = []
        self.vector = vector

    def add_operator(self, op):
        self.operators.append(op)

    def add_regulator(self, reg):
        self.regulators.append(reg)

    def add_reporter(self, rep):
        self.reporters.append(rep)

    def step(self, growth_rate=1, t=0, dt=0.1):
        for op in self.operators:
            expression_rate = op.expression_rate(t, dt)
            op.output.express(expression_rate)

        for regulator in self.regulators:
            regulator.step(growth_rate, dt)

        for reporter in self.reporters:
            reporter.step(growth_rate, dt)


        
