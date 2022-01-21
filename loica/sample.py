class Sample:
    """
    Representation of a sample that encapsulates GeneticNetwork and Metabolism.
    Incorporate environment information such as Supplements or chemicals, strain and media. 
    Ex: 1 well in a plate, single cell.
    ...

    Attributes
    ----------
    genetic_network : GeneticNetwork
        genetic network that is part of the sample
    metabolism : Metabolism
        metabolism that drives the genetic network in the sample
    assay : Assay
        assay to which this sample belongs
    media : str
        Name of the media in the sample
    strain : str
        Name of the strain in the sample
    
     Methods
    -------
    add_supplement(supplement, concentration)
        stablishes the concentration of Supplement
    """
    def __init__(self, 
            genetic_network=None, 
            metabolism=None, 
            assay=None,
            media=None,
            strain=None,
            ):
        self.genetic_network = genetic_network
        self.metabolism = metabolism
        self.media = media
        self.strain = strain
        self.vector = self.genetic_network.vector
        if self.genetic_network:
            self.reporters = self.genetic_network.reporters
        if metabolism:
            self.biomass = self.metabolism.biomass
        self.supplements = {}

    def initialize(self):
        self.genetic_network.initialize()

    def set_supplement(self, supplement, concentration):
        self.supplements[supplement] = concentration

    def set_regulator(self, name, concentration):
        for reg in self.genetic_network.regulators:
            if reg.name == name:
                reg.init_concentration = concentration
            else: pass

    def set_reporter(self, name, concentration):
        for rep in self.genetic_network.reporters:
            if rep.name == name:
                rep.init_concentration = concentration
            else: pass

    def step(self, t, dt, stochastic=False):
        if self.genetic_network and self.metabolism:
            growth_rate = self.metabolism.growth_rate(t)
            for supp,conc in self.supplements.items():
                supp.concentration = conc
            if stochastic:
                self.genetic_network.step_stochastic(growth_rate, t, dt)
            else:
                self.genetic_network.step(growth_rate, t, dt)
            self.reporters = self.genetic_network.reporters


