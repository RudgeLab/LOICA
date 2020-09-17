class Sample:
    def __init__(self, 
            circuit=None, 
            metabolism=None, 
            init_biomass=0,
            assay=None,
            media=None,
            strain=None,
            vector=None
            ):
        self.circuit = circuit
        self.metabolism = metabolism
        self.biomass = init_biomass
        self.media = media
        self.strain = strain
        self.vector = vector
        if self.circuit:
            self.reporters = self.circuit.reporters

    def step(self, t, dt):
        if self.circuit and self.metabolism:
            growth_rate = self.metabolism.growth_rate(t)
            profile = self.metabolism.profile(t)
            self.circuit.step(growth_rate, dt)
            self.reporters = self.circuit.reporters
        self.biomass = self.metabolism.biomass(t)


