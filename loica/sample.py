class Sample:
    def __init__(self, circuit=None, metabolism=None, init_biomass=0):
        self.circuit = circuit
        self.metabolism = metabolism
        self.signals = {reg.name:reg.concentration for reg in self.circuit.regulators}
        self.biomass = init_biomass

    def step(self, t, dt):
        if self.circuit and self.metabolism:
            growth_rate = self.metabolism.growth_rate(t)
            profile = self.metabolism.profile(t)
            self.circuit.step(growth_rate, dt)
            self.signals = {reg.name:reg.concentration for reg in self.circuit.regulators}
        self.biomass = self.metabolism.biomass(t)


