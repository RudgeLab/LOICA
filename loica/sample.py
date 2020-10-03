class Sample:
    def __init__(self, 
            circuit=None, 
            metabolism=None, 
            assay=None,
            media=None,
            strain=None,
            ):
        self.circuit = circuit
        self.metabolism = metabolism
        self.media = media
        self.strain = strain
        self.vector = self.circuit.vector
        if self.circuit:
            self.reporters = self.circuit.reporters
        if metabolism:
            self.biomass = self.metabolism.biomass
        self.supplements = {}

    def initialize(self):
        self.circuit.initialize()

    def add_supplement(self, supplement, concentration):
        self.supplements[supplement] = concentration

    def step(self, t, dt):
        if self.circuit and self.metabolism:
            growth_rate = self.metabolism.growth_rate(t)
            for supp,conc in self.supplements.items():
                supp.concentration = conc
            self.circuit.step(growth_rate, t, dt)
            self.reporters = self.circuit.reporters


