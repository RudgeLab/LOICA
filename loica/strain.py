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
    cell_volume : int | float
        average cell volume. 
        Default cell volume has been taken from Supplementary data in:
            Yáñez Feliú, G., Vidal, G., Muñoz Silva M., Rudge, T. J. (2020). 
            Novel Tunable Spatio-Temporal Patterns From a Simple Genetic Oscillator Circuit. 
            Frontiers in Bioengineering and Biotechnology 8, 893. 
            https://doi.org/10.3389/FBIOE.2020.00893/BIBTEX
        
    strain : str
        Name of the strain in the sample (link to Flapjack)
    
    """

    def __init__(self, name=None, genetic_network=None, metabolism=None, cell_volume=10**-15, strain=None):
        self.name = name
        self.genetic_network = genetic_network
        self.metabolism = metabolism
        self.cell_volume = cell_volume # default volume is in L
        self.cell_number = 0
        if self.genetic_network:
            self.reporters = self.genetic_network.reporters
            self.regulators = self.genetic_network.regulators
            self.gene_products = self.genetic_network.gene_products
        if metabolism:
            self.biomass = self.metabolism.biomass
            self.growth_rate = self.metabolism.growth_rate

    # TODO: add functions add_reporter and etc
    
        

        