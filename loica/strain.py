class Strain:
    """
    Representation of a strain that encapsulates GeneticNetwork and Metabolism.
    ...

    Attributes
    ----------
    name : str
        name of the strain
    genetic_network : GeneticNetwork
        genetic network that is part of the strain
    metabolism : Metabolism
        metabolism that drives the genetic network in the strain
    strain : str
        Name of the strain in the sample (link to Flapjack)
    
    """

    def __init__(self, name=None, genetic_network=None, metabolism=None, strain=None):
        self.name = name
        self.genetic_network = genetic_network
        self.metabolism = metabolism
        if self.genetic_network:
            self.reporters = self.genetic_network.reporters
            self.gene_products = self.genetic_network.regulators + self.reporters
        if metabolism:
            self.biomass = self.metabolism.biomass
            self.growth_rate = self.metabolism.growth_rate

    # TODO: add functions add_reporter and etc
    
        

        