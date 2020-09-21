class Sample:
    def __init__(self, 
            circuit=None, 
            metabolism=None, 
            init_biomass=0,
            assay=None,
            media=None,
            strain=None,
            ):
        self.circuit = circuit
        self.metabolism = metabolism
        self.biomass = init_biomass
        self.media = media
        self.strain = strain
        self.vector = self.circuit.vector
        if self.circuit:
            self.reporters = self.circuit.reporters
        self.supplements = {}

    def add_supplement(self, supplement, concentration):
        self.supplements[supplement] = concentration

    def step(self, t, dt):
        if self.circuit and self.metabolism:
            growth_rate = self.metabolism.growth_rate(t)
            for supp,conc in self.supplements.items():
                supp.concentration = conc
            self.circuit.step(growth_rate, t, dt)
            self.reporters = self.circuit.reporters
        self.biomass = self.metabolism.biomass(t)


