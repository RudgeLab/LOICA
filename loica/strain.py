from .sample import Sample


class Strain(Sample):
    """
    Representation of a sample that contains only one strain.
    Encapsulates GeneticNetwork and Metabolism. Incorporate environment information 
    such as Supplements or chemicals, strain and media. 

    Child of Sample
    
    """

    def __init__(self, 
            genetic_network=None, 
            metabolism=None, 
            assay=None,
            media=None,
            strain=None,
            ):
        super().__init__(genetic_network, metabolism, assay, media, strain)  
        # Do I need code below?  
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


