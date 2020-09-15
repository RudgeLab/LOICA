class Sample:
    def __init__(self, circuit=None, metabolism=None, init_biomass=0):
        self.circuit = circuit
        self.metabolism = metabolism
        self.signals = {n.name: n.concentration for n in self.circuit.nodes}
        self.biomass = init_biomass

    def step(self, t, dt):
        if self.circuit and self.metabolism:
            growth_rate = self.metabolism.growth_rate(t)
            profile = self.metabolism.profile(t)
            self.circuit.step(growth_rate, dt)
            self.signals = {n.name: n.concentration for n in self.circuit.nodes}
        self.biomass = self.metabolism.biomass(t)


