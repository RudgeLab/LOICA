class Sample:
    """
    Representation of a sample that encapsulates GeneticNetwork and Metabolism.
    Incorporate environment information such as Supplements or chemicals, strain and media. 
    Ex: 1 well in a plate, single cell.
    ...

    Attributes
    ----------
    genetic_network : List[GeneticNetwork] or GeneticNetwork
        genetic network that is part of the sample
    metabolism : list[Metabolism] or Metabolism
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
        # TODO: fix
        # self.vector = self.genetic_network.vector
        self.reporters = []
        self.gene_products = []

        # adding all reporters into list
        if type(self.genetic_network)==list:
            for genetic_network in self.genetic_network:
                self.reporters += genetic_network.reporters
                self.gene_products += genetic_network.reporters + genetic_network.regulators
        else:
            self.reporters = self.genetic_network.reporters
            self.gene_products = self.reporters + self.genetic_network.regulators
        # sort by name so gene products with the same identity would be together
        self.gene_products.sort(key=lambda x: x.name)

        # TODO: add code that throws an error if extracelluar degradation rates of 
        # geneproducts with the same identity are different

        if self.metabolism:
            if type(self.metabolism)==list:
                self.biomass = []
                for m in metabolism:
                    self.biomass.append(m.biomass)
            else:       
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
                if type(self.genetic_network)==list and type(self.metabolism)==list:
                    for (gn, b, g_rate) in zip(self.genetic_network, biomass, growth_rate):
                        gn.step(b, g_rate, t, dt)
                elif type(self.genetic_network)==list:
                    for gn in self.genetic_network:
                        gn.step(biomass, growth_rate, t, dt) 
                else:
                    self.genetic_network.step(biomass, growth_rate, t, dt)
                
                # add degradation of extracellular signal
                # TODO: but here is an issue. The extracellular concentration should 
                # update once it changes, therefore during each iteration of 
                # genetic_network.step. There should be something like 
                # for gp in gn.gene_products:
                #   if gp.name in self.gene_products.name:
                #       all gene products with the same name get the concentration 
                #       updated.
                # TODO: make this a function and then stick after each step?
                for i, gp in enumerate(self.gene_products):
                    if gp in gn.gene_product:
                        if gp.name == self.gene_products[i-1]:
                            self.gene_products[i-1].ext_conc = gp.ext_conc
                        if gp.name == self.gene_products[i+1]:
                            self.gene_products[i+1].ext_conc = gp.ext_conc
                #  
                for i, gene_product in enumerate(self.gene_products):
                    if gene_product.name != self.gene_products[i-1].name:
                        ext_degr = gene_product.ext_conc * gene_product.ext_degr_rate
                        new_ext_conc = gene_product.ext_conc - ext_degr
                        gene_product.ext_conc = new_ext_conc
                    else: 
                        gene_product.ext_conc = self.gene_products[i-1].ext_conc
                        self.gene_products[i-1].ext_conc = gene_product.ext_conc

                # option for changing external concentration based on the gp.ext_difference
                for i, gp in enumerate(self.gene_products):
                    # problem - will apply degradation rate at each step
                    ext_degr = gp.ext_conc * gp.ext_degr_rate
                    new_ext_conc = gp.ext_conc - ext_degr + gp.ext_difference
                    gp.ext_conc = new_ext_conc
                    if gp.name == self.gene_products[i-1].name:
                        new_ext_conc = gp.ext_conc + gp.ext_difference
                        gp.ext_conc = new_ext_conc
                        self.gene_products[i-1].ext_conc = new_ext_conc
                    if gp.name == self.gene_products[i+1].name:
                        new_ext_conc = gp.ext_conc + gp.ext_difference
                        gp.ext_conc = new_ext_conc
                        self.gene_products[i+1].ext_conc = new_ext_conc  
                    #  problem - will update only adjacent two 

                
                    

            # if type(self.genetic_network)==list:
            #     for gn in self.genetic_network:
            #         for reporter in gn.reporters:
            #             self.reporters.append(reporter)
            # else:
            #     self.reporters = self.genetic_network.reporters


