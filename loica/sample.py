from itertools import groupby
import numpy as np
from random import shuffle

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
                growth_rate = []
                biomass = []
                for m in self.metabolism:
                    growth_rate.append(m.growth_rate)
                    biomass.append(m.biomass)
            else:
                growth_rate = self.metabolism.growth_rate
                biomass = self.metabolism.biomass

        self.supplements = {}

        # create an options list for stochastic simulation 
        # (used in self.total_substep_stochastic(self))
        if type(self.genetic_network)==list:
            self.options = self.genetic_network 
        else:
            self.options = []
            self.options.append(self.genetic_network)
        self.options.append("extracellular space")

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
    # TODO: work on it
    # def supplement_is_gp(self, supplement):
    #     for gp in self.extracellular_space:
    #         if supplement.name == gp.name:
    #             gp.extracellular_conc=supplement.concentration

    # TODO: update these two methods
    # def set_regulator(self, name, concentration):
    #     for reg in self.genetic_network.regulators:
    #         if reg.name == name:
    #             reg.init_concentration = concentration
    #         else: pass

    # def set_reporter(self, name, concentration):
    #     for rep in self.genetic_network.reporters:
    #         if rep.name == name:
    #             rep.init_concentration = concentration
    #         else: pass

    def external_step(self):
        """ 
            method that calculates the change in the extracellular concentration
            due to degradation

            deterministic
        """
        if type(self.gene_products[0])==list:
            for group in self.gene_products: 
                ext_degr = group[0].ext_conc * group[0].ext_degr_rate
                new_ext_conc = group[0].ext_conc - ext_degr
                for gp in group:
                    gp.ext_conc = new_ext_conc
                # test
                # print(f'New ext_conc of {group[0].name} = {new_ext_conc}')
        else:
            for gp in self.gene_products:
                ext_degr = gp.ext_conc * gp.ext_degr_rate
                new_ext_conc = gp.ext_conc - ext_degr
                gp.ext_conc = new_ext_conc

    def update_ext_conc(self):
        '''
            Method to update external concentration of all gene products based on the
            geneproduct.ext_difference 
        '''
        if type(self.gene_products[0])==list:
            for group in self.gene_products:
                concentration_change = 0
                for gp in group:
                    concentration_change += gp.ext_difference   
                new_ext_conc = group[0].ext_conc + concentration_change
                for gp in group:
                    gp.ext_conc = new_ext_conc
        else:
            for gp in self.gene_products:
                new_ext_conc = gp.ext_conc + gp.ext_difference
                gp.ext_conc = new_ext_conc

    def st_external_substep(self, tau_=None):
        """ 
            method similar to GeneticNetwork.substep()
            but only for extracellular degradation

            stochastic
        """
        # Propensities
        a = []

        # Compute propensities for degradation of gene products
        if type(self.gene_products[0])==list:
            for group in self.gene_products:
                # degradation reaction
                a.append(group[0].ext_degr_rate * group[0].ext_conc)
        else:
            for gp in self.gene_products:
                # degradation reaction
                a.append(gp.ext_degr_rate * gp.ext_conc)
        
        # Make list of propensities into array
        a = np.array(a)
        # Total of propensities
        A = a.sum()
        
        # Time step
        if tau_:
            tau = tau_
        else:
            tau = 1/A * np.log(1/np.random.random())

        # Random number to select next reaction
        a_i = np.random.random() * A

        # Find reaction and update gene product levels
        for i, group in enumerate(self.gene_products):
            if a_i < np.sum(a[:i+1]):
                # Extracellular degradation of geneproduct gp
                for gp in group:
                    gp.ext_conc -= 1
                break

        # Return elapsed time if it was not predefined
        if not tau_:
            return tau
        
    def correct_metabolism(self, genetic_network):
        ''' method to pick correct growth rate and biomass '''
        if type(self.metabolism)==list:
            i = self.genetic_network.index(genetic_network)
            biomass = self.biomass[i]
            growth_rate = self.growth_rate[i]
        else:
            biomass = self.biomass
            growth_rate = self.growth_rate       
        return growth_rate, biomass
        
    def total_substep_stochastic(self, t=0, dt=0.1):
        '''
            method that links stochastic substeps for each genetic network and 
            extracellular space
        '''
        # shuffle the listbiomass, t, dt, growth_rate
        shuffle(self.options)

        # get tau by running substep for the first item in the shuffled list
        # and then use this tau in other substeps
        if self.options[0]=="extracellular space":
            tau = self.st_external_substep()
            for gn in self.options[1:len(self.options)]:
                growth_rate = self.correct_metabolism(gn)[0]
                biomass = self.correct_metabolism(gn)[1]                     
                gn.substep_stochastic(t, dt, growth_rate, biomass, tau)
            self.update_ext_conc()
        else:
            growth_rate = self.correct_metabolism(self.options[0])[0]
            biomass = self.correct_metabolism(self.options[0])[1]
            tau = self.options[0].substep_stochastic(t, dt, growth_rate, biomass)
            for gn in self.options[1:len(self.options)]:
                if gn == "extracellular space":
                    self.st_external_substep(tau_=tau)
                    self.update_ext_conc()
                else:
                    growth_rate = self.correct_metabolism(gn)[0]
                    biomass = self.correct_metabolism(gn)[1] 
                    gn.substep_stochastic(t, dt, growth_rate, biomass, tau)
                    self.update_ext_conc()

    def step_stochastic(self, t=0, dt=0.1):
        ''' 
            same as GeneticNetwork.step_stochastic(), but uses either 
            self.total_substep_stochastic() or GeneticNetwork.step_stochastic()
        '''
        delta_t = 0
        if type(self.genetic_network)==list:
            # if many gene networks - use total_substep_stochastic
            while delta_t < dt:
                #print(f'Elapsed time: {delta_t}')
                delta_t += self.total_substep_stochastic(t, dt)
        else:
            for gp in self.gene_products:
                if gp.diffusion_rate != 0:
                    while delta_t < dt:
                        #print(f'Elapsed time: {delta_t}')
                        delta_t += self.total_substep_stochastic(t, dt)
                    break
            if delta_t == 0:
                self.genetic_network.step_stochastic(t, dt, self.growth_rate)

    def step(self, t, dt, stochastic=False):
        if self.genetic_network and self.metabolism:
            for supp,conc in self.supplements.items():
                supp.concentration = conc

            if stochastic:
                self.step_stochastic(t, dt)
            else:
                if type(self.genetic_network)==list and type(self.metabolism)==list:
                    for (gn, b, g_rate) in zip(self.genetic_network, self.biomass, self.growth_rate):
                        # test
                        print(f'In network {gn} at biomass {b}')
                        gn.step(b, g_rate, t, dt)
                elif type(self.genetic_network)==list:
                    for gn in self.genetic_network:
                        gn.step(self.biomass, self.growth_rate, t, dt) 
                else:
                    self.genetic_network.step(self.biomass, self.growth_rate, t, dt)
                # update the exctracellular concentration
                self.external_step()
                self.update_ext_conc()



