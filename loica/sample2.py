class Sample:
    """
    Representation of a sample that encapsulates GeneticNetwork and Metabolism.
    Incorporate environment information such as Supplements or chemicals, strain and media. 
    Ex: 1 well in a plate, single cell.
    ...

    Attributes
    ----------
    genetic_network : List[GeneticNetwork]
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
    # TODO: determine whether I need different metabolism for each strain (in general)
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
        self.reporters = []

        # adding all reporters into list
        if type(self.genetic_network)==list:
            for genetic_network in self.genetic_network:
                for reporter in genetic_network.reporters:
                    self.reporters.append(reporter)
        else:
            self.reporters = self.genetic_network.reporters

        
        # setting up extracellular space
        # self.extr_conc = []
        # for gn in self.genetic_network:
        #     genetic_producs = gn.reporters + gn.regulators

            """ 
            make an array or dictionary which looks like this:

            {gp:conc, gp2:conc2} - dictionary. 
            [conc1, conc2, conc3, etc]
            
            however, I need to make sure that:
            a) gene products with the same identity are not repeated
            I could add attribute "extracellular_conc" to GeneProduct, set it to 0 by 
            default. Then when gene products from different strains are "merged", 

            maybe:
            [ [[gene products with the same identity], extracellular concentration],
              [[gene products with the same identity], extracellular concentration],
                ]  

            at the moment I am not sure what is the best way to go about this, I need
            to decide later.
            
            """

        # this assumes that metabolism is the same for all 3 strains.
        # might want to change it in the future, and code interactions between strains
        # into metabolism.py
        if metabolism:
            self.biomass = self.metabolism.biomass

        self.supplements = {}

    def initialize(self):
        if type(self.genetic_network)==list:
            for gn in self.genetic_network:
                gn.initialize()
        else:
            self.genetic_network.initialize()
        # and initialise extracellular concentrations if needed

    def set_supplement(self, supplement, concentration):
        self.supplements[supplement] = concentration

    # add function that checks whether supplement is the same as gene product
    # if yes, add the supplement concentration to the concentration of gp in the
    # extracellular space
    def supplement_is_gp(self, supplement):
        for gp in self.extracellular_space:
            if supplement.name == gp.name:
                gp.extracellular_conc=supplement.concentration

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
            # check if metabolism is list
            if type(self.metabolism)==list:
                growth_rate = []
                biomass = []
                for m in self.metabolism:
                    growth_rate.append(m.growth_rate(t))
                    biomass.append(m.biomass(t))
            else:
                growth_rate = self.metabolism.growth_rate(t)
                biomass = self.metabolism.biomass(t)

            for supp,conc in self.supplements.items():
                supp.concentration = conc

            if stochastic:
                # I need for all strains to have steps either simultaneously or
                # randomly, with the simultaneous change extracellularly
                if type(self.genetic_network)==list:
                    for gn in self.genetic_network:
                        gn.step_stochastic(growth_rate, t, dt)
                else:
                    self.genetic_network.step_stochastic(growth_rate, t, dt)
            else:
                if type(self.genetic_network)==list and type(biomass)==list:
                    for (gn, b) in zip(self.genetic_network, biomass):
                        gn.step(b, growth_rate, t, dt)
                elif type(self.genetic_network)==list:
                    for gn in self.genetic_network:
                        gn.step(biomass, growth_rate, t, dt) 
                else:
                    self.genetic_network.step(biomass, growth_rate, t, dt)
            
            if type(self.genetic_network)==list:
                for gn in self.genetic_network:
                    for reporter in gn.reporters:
                        self.reporters.append(reporter)
            else:
                self.reporters = self.genetic_network.reporters


