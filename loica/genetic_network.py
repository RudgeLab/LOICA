class GeneticNetwork():
    def __init__(self):
        self.tus = []
        self.regulators = []

    def add_tu(self, tu):
        self.tus.append(tu)

    def add_regulator(self, reg):
        self.regulators.append(reg)

    def step(self, growth_rate=1, dt=0.1):
        for tu in self.tus:
            expression_rate = tu.expression_rate()
            tu.output.express(expression_rate)

        for regulator in self.regulators:
            regulator.step(growth_rate, dt)



        
