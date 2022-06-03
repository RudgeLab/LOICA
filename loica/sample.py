from itertools import groupby

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
        
        '''
        add resources - which are depleted by all cells, thus limiting growth.
        it could be something like

        parameter (resources = 1000000)
        self.resources = resources

        # this can be put in metabolism
        Metabolism(consumption_rate=2)
        def deplete(self, resources)
            consumption = self.biomass * self.consumption_rate
            new_resources = resources - consumption
            resources = new_resources
            if resources <= 0:
                self.growth_rate = 0
        
        so each strain can have different consumption rate

        and then it can be incorporated here or in assay.

        TODO: think about it further. most likely will be in assay
        '''

        # adding all reporters into list
        if type(self.genetic_network)==list:
            list_gp = []
            for genetic_network in self.genetic_network:
                self.reporters += genetic_network.reporters
                list_gp += genetic_network.reporters + genetic_network.regulators
            # sort by name so gene products with the same identity would be together
            list_gp.sort(key=lambda x: x.name)
            # group into sublists based on the identity and add to self.gene_products
            for name, identical_gproducts in groupby(list_gp, lambda x: x.name):
                self.gene_products.append(list(identical_gproducts))
        else:
            self.reporters = self.genetic_network.reporters
            self.gene_products = self.reporters + self.genetic_network.regulators

        if self.metabolism:
            if type(self.metabolism)==list:
                self.biomass = []
                for m in metabolism:
                    self.biomass.append(m.biomass)
            else:       
                self.biomass = self.metabolism.biomass

        self.supplements = {}

    def set_extracel_degr(self, chemical_name, ext_degr_rate):
        ''' 
            this method ensures that all gene products with the same identity 
            have the same extracellular degradation rate
        '''
        if type(self.gene_products[0])==list:
            for group in self.gene_products:
                if group[0].name == chemical_name:
                    for gp in group:
                        gp.ext_degr_rate = ext_degr_rate
        else:
            for gp in self.gene_products:
                if gp.name == chemical_name:
                    gp.ext_degr_rate = ext_degr_rate

    def initialize(self):
        if type(self.genetic_network)==list:
            for gn in self.genetic_network:
                gn.initialize()
        else:
            self.genetic_network.initialize()

    def set_supplement(self, supplement, concentration):
        self.supplements[supplement] = concentration

    # add function that checks whether supplement is the same as gene product
    # if yes, add the supplement concentration to the concentration of gp in the
    # extracellular space
    def supplement_is_gp(self, supplement):
        for gp in self.extracellular_space:
            if supplement.name == gp.name:
                gp.extracellular_conc=supplement.concentration

    # TODO: update these two methods
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

    def external_step(self):
        """ 
        method that calculates the total change in the extracellular concentration
        and updates it

        deterministic
        """
        if type(self.gene_products[0])==list:
            for group in self.gene_products:
                concentration_change = 0
                for gp in group:
                    concentration_change += gp.ext_difference   
                ext_degr = group[0].ext_conc * group[0].ext_degr_rate
                new_ext_conc = group[0].ext_conc - ext_degr + concentration_change
                for gp in group:
                    gp.ext_conc = new_ext_conc
                # test
                print(f'New ext_conc of {group[0].name} = {new_ext_conc}')
        else:
            for gp in self.gene_products:
                ext_degr = gp.ext_conc * gp.ext_degr_rate
                new_ext_conc = gp.ext_conc - ext_degr + gp.ext_difference
                gp.ext_conc = new_ext_conc


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
                        # test
                        print(f'In network {gn} at biomass {b}')
                        gn.step(b, g_rate, t, dt)
                elif type(self.genetic_network)==list:
                    for gn in self.genetic_network:
                        gn.step(biomass, growth_rate, t, dt) 
                else:
                    self.genetic_network.step(biomass, growth_rate, t, dt)
                # update the exctracellular concentration
                self.external_step()



