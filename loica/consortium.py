from .sample import Sample


class Consortium(Sample):
    """
    Representation of microbial consortium (a child of Sample)
    Consortium contains GeneticNetworks (which equal strains) and their respective 
    Metabolism. The strains are connected through extracellular space/medium.
    
    Attributes
    ----------
    strains : List[GeneticNetwork]
        strains that form consortium
    metabolism : List[Metabolism]
        metabolism that drives the genetic network of each strain
    media : str
        Name of the media in the consortium
    strain : List[str]
        Name of the strains in the consortium
    
    """

    """ 
    I need to have something similar to metabolism here
    So there will be either nested list or dictionary containing reporters and their 
    extracelllar concentration, or I need to make gene products to have extracellular
    concentration

    For example, extracellular_conc = [[reporter1, conc1], [reporter2, conc2]]
    And each extracellular concentration will be set to 0 in the beginning
    extracellular_conc = [[reporter1, 0], [reporter2, 0]]

    Since supplement can be a reporter, I need to change reporter concentration to be 
    equal to supplement

    so let's say supplement = lc.Supplement(name="Reporter1", concentration=5)
    and reporter1.name = "Reporter1"
    since supplement.name == reporter1.name
    extracellular_conc = [[reporter1, 5], [reporter2, 0] at time 0
    """
    def __init__(self,
            strains=None, 
            metabolism=None, 
            media=None,
            strain=None,
            ):
        super().__init__(metabolism, media)
        self.strains = strains
        self.strain = strain
        #self.vector = self.genetic_network.vector

        if self.strains:
            # code that pulls reporters from all strains, and removes repeating
            # add code for latter (not sure if I need to remove repeating)
            self.reporters = []
            for strain in self.strains:
                self.reporters.append(strain.reporters)
            # removing duplicates
            self.reporters = list(dict.fromkeys(self.reporters))
        if metabolism:
            self.biomass = []
            for strain_metabolism in self.metabolism:
                self.biomass.append(strain_metabolism.biomass)
        self.supplements = {} 

        self.extracellular_conc = {}


    # function that checks whether supplement is the same as reporter
    # if yes, set the extracellular concentration of the reporter to be equal 
    # to supplemet
    def supplement_is_reporter(self, supplement):
        for reporter in self.reporters:
            if supplement.name == reporter.name:
                reporter.extracellular_conc=supplement.concentration

    # method that sets all concentrations to zero both within and outside cells
    def initialize_both(self):
        for strain in self.strains:
            strain.initialize()
            for reporter in strain.reporters
