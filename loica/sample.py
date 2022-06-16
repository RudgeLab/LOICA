from itertools import groupby
import numpy as np
from random import shuffle
from .strain import Strain

class Sample:
    """
    Representation of a sample that encapsulates either one Strain or multiple (consortium).
    Incorporate environment information such as Supplements or chemicals and media. 
    Ex: 1 well in a plate, single cell.
    ...

    Attributes
    ----------
    strain : List[Strain] or Strain
        strain that is part of the sample
    media : str
        Name of the media in the sample
    
     Methods
    -------
    add_supplement(supplement, concentration)
        establishes the concentration of Supplement
    """
    def __init__(self, strain=None, media=None): # resources=None

        self.strain = []
        self.media = media
        self.reporters = []
        self.gene_products = []
        self.growth_rate = []
        self.biomass = []
        self.supplements = {}
        # self.resources = resources        

        if issubclass(type(strain), Strain):
            self.strain = [strain]
        elif type(strain)==list:
            for s in strain:
                if issubclass(type(s), Strain):
                    self.strain.append(s)
                else: print('Unsupported Type, it should be a Strain')
        
        # adding all reporters, gene products, growth rates and biomass to respective 
        # lists
        if self.strain:
            list_gp = []
            for s in self.strain:
                if s.genetic_network:
                    self.reporters += s.reporters
                    list_gp += s.gene_products
                # assign strain to each gene_product (for self.catch_negative_conc())
                for gp in s.gene_products:
                    gp.strain = s
                if s.metabolism:
                    self.growth_rate.append(s.growth_rate)
                    self.biomass.append(s.biomass)
            # sort by name so gene products with the same identity would be together
            list_gp.sort(key=lambda x: x.name)
            # group into sublists based on the identity and add to self.gene_products
            for name, identical_gproducts in groupby(list_gp, lambda x: x.name):
                self.gene_products.append(list(identical_gproducts))

        '''
        # TODO: adapt for using Strain
        # create an options list for stochastic simulation 
        # (used in self.total_substep_stochastic(self))
        if type(self.genetic_network)==list:
            # copy list without further changing original
            # slicing is fastest method according to @cryo at 
            # https://stackoverflow.com/questions/2612802/how-do-i-clone-a-list-so-that-it-doesnt-change-unexpectedly-after-assignment
            self.options = self.genetic_network[:]
        else:
            self.options = []
            self.options.append(self.genetic_network)

        '''

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
        
        # this ensures that self.st_external_substep() is called later in self.total_substep_stochastic()
        if "extracellular space" not in self.options:
            self.options.append("extracellular space")

    def initialize(self):
        for s in self.strain:
            s.genetic_network.initialize()

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

    def external_step(self, dt):
        """ 
            method that calculates the change in the extracellular concentration
            due to degradation

            deterministic
        """
        for group in self.gene_products: 
            ext_degr = group[0].ext_conc * group[0].ext_degr_rate
            new_ext_conc = group[0].ext_conc - ext_degr * dt
            for gp in group:
                gp.ext_conc = new_ext_conc

    def update_ext_conc(self, stochastic=False):  
        '''
            Method to update external concentration of all gene products based on the
            geneproduct.ext_difference 
        '''
        for group in self.gene_products:
            concentration_change = 0
            for gp in group:
                concentration_change += gp.ext_difference   
            new_ext_conc = group[0].ext_conc + concentration_change - group[0].ext_degraded
            if new_ext_conc<0:
                new_ext_conc = self.catch_negative_conc(group, stochastic)
            for gp in group:
                gp.ext_conc = new_ext_conc

    def catch_negative_conc(self, group, stochastic=False):
        ''' 
            used if external concentration becomes negative due to multiple strains
            taking molecules out of the extracellular space simultaneously.
            randomly updates extracellular concentration
            if result is negative, diffusion does not happen and molecules
            are returned to the cell
        '''
        if stochastic:
            ext_options = ["degradation"]
            for gp in group:
                ext_options.append(gp) 
            shuffle(ext_options)
            new_ext_conc = group[0].ext_conc
            for opt in ext_options:
                if opt == "degradation":
                    new_ext_conc -= group[0].ext_degraded
                    if new_ext_conc < 0:
                        new_ext_conc = 0
                else:
                    # this chunk is not working as it should
                    new_ext_conc += opt.ext_difference
                    if new_ext_conc < 0:
                        new_ext_conc -= opt.ext_difference
                        opt.concentration -= 1 # causes opt.concentration to become negative for some reason
                        # test
                        # if opt.concentration<0:
                        #     print(f"Int conc of {opt.name} {opt} is {opt.concentration}")
            return new_ext_conc
        else:
            # list of gene products that diffuse out of the extracellular space
            diffused_out = []
            for gp in group:
                if gp.ext_difference<0:
                    diffused_out.append[gp]
                else:
                    new_ext_conc += gp.ext_difference
            ideal_diffusion_out = 0
            for gp in diffused_out:
                # calculate what would be the change of concentration if there were 
                # enough resources
                ideal_diffusion_out -= gp.ext_difference
            for gp in diffused_out:
                # calculate proportion of molecules each strain would take in from total
                # then multiply by the available concentration to get the concentration
                # change that happened
                can_take = (ideal_diffusion_out/(-gp.ext_difference)*new_ext_conc)
                # calculate "extra" concentration of the gp in the strain it belongs to
                extra_taken = (-gp.ext_difference-can_take)/gp.strain.biomass
                # correct the internal gp concentration
                fixed_conc = gp.concentration - extra_taken
                gp.concentration = fixed_conc
            return 0
                
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
                # reset how much degraded
                group[0].ext_degraded = 0
        else:
            for gp in self.gene_products:
                # degradation reaction
                a.append(gp.ext_degr_rate * gp.ext_conc)
                # reset how much degraded
                gp.ext_degraded = 0
        
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
                # Mark extracellular degradation of geneproduct in group
                group[0].ext_degraded = 1
                # #test
                # print(f'{group[0].name} ext degradation')
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
            for gn in self.options[1:]:
                growth_rate = self.correct_metabolism(gn)[0](t)
                biomass = self.correct_metabolism(gn)[1](t)                    
                gn.substep_stochastic(t, dt, growth_rate, biomass, tau)
            self.update_ext_conc()
        else:
            growth_rate = self.correct_metabolism(self.options[0])[0](t)
            biomass = self.correct_metabolism(self.options[0])[1](t)
            tau = self.options[0].substep_stochastic(t, dt, growth_rate, biomass)
            for gn in self.options[1:]:
                if gn == "extracellular space":
                    self.st_external_substep(tau_=tau)
                else:
                    growth_rate = self.correct_metabolism(gn)[0](t)
                    biomass = self.correct_metabolism(gn)[1](t)
                    gn.substep_stochastic(t, dt, growth_rate, biomass, tau)
            self.update_ext_conc()
                    
        return tau

    def substep_stochastic(self, t=0, dt=0.1, growth_rate=1, biomass=1):
        ''' fully stochasic, only one reaction happens out of all'''
        # Propensities
        a = []

        # Compute expression rates
        if type(self.genetic_network)==list:
            for gn in self.genetic_network:
                for op in gn.operators:
                    expression_rate = op.expression_rate(t, dt)
                    output = op.output
                    if type(op.output)!=list:
                        output = [output]
                    for o in output:
                        o.express(expression_rate)
        # TODO: add code for non-list

        # Compute propensities for production, degradation, diffusion in and out of gene products
        if type(self.gene_products[0])==list:
            for group in self.gene_products:
                for gp in group:
                    # Production reaction
                    a.append(gp.expression_rate)
                    # Degradation
                    a.append((gp.degradation_rate + growth_rate) * gp.concentration)
                    # Diffusion out of cell
                    a.append(gp.diffusion_rate*gp.concentration)
                    # Difusion into the cell
                    a.append(gp.diffusion_rate*gp.ext_conc)
                    # External degradation
                    a.append(gp.ext_conc*gp.ext_degr_rate)
                    # reset values
                    gp.ext_difference = 0
            
        # TODO: add other instances

        # Make list of propensities into array
        a = np.array(a)
        # Total of propensities
        A = a.sum()
        
        # Time step
        tau = 1/A * np.log(1/np.random.random())
        # Random number to select next reaction
        a_i = np.random.random() * A

        # Find reaction and update gene product levels
        if type(self.gene_products[0])==list:
            x = 0
            complete = False
            for group in self.gene_products:
                x += len(group)
                for ii, gp in enumerate(group):
                    i = ii+x
                    if a_i < np.sum(a[:i*5+1]):
                        # Production of geneproduct gp
                        gp.concentration += 1
                        complete = True
                        break
                    elif a_i < np.sum(a[:i*5+2]):
                        # Degradation of geneproduct gp
                        gp.concentration -= 1
                        complete = True
                        break
                    elif a_i < np.sum(a[:i*5+3]):
                        # Diffusion of geneproduct gp out of cell
                        gp.concentration -= 1
                        gp.ext_difference = biomass
                        complete = True
                        break
                    elif a_i < np.sum(a[:i*5+4]):
                        # Diffusion of geneproduct gp into the cell
                        gp.concentration += 1
                        gp.ext_difference = - biomass
                        complete = True
                        break
                    elif a_i < np.sum(a[:i*5+5]):
                        # External degradation of geneproduct gp
                        gp.ext_difference -= 1
                        complete = True
                        break
                if complete:
                    break
        
        # Reset expression rates for next step and update external concentration
        if type(self.gene_products[0])==list:
            for group in self.gene_products:
                concentration_change = 0
                for gp in group:
                    gp.expression_rate = 0
                    concentration_change += gp.ext_difference   
                new_ext_conc = group[0].ext_conc + concentration_change
                for gp in group:
                    gp.ext_conc = new_ext_conc
         # TODO: add other cases
        
        # Return elapsed time
        return tau

    def step_stochastic(self, t=0, dt=0.1, growth_rate=1, biomass=1):
        ''' 
            same as GeneticNetwork.step_stochastic(), but uses either 
            self.total_substep_stochastic() or GeneticNetwork.step_stochastic()
        '''
        delta_t = 0
        while delta_t < dt:
                # print(f'Elapsed time: {delta_t}')
                delta_t += self.substep_stochastic(t, dt, growth_rate, biomass)
        ''' 
        # old code for semi-stochastic and fully stochastic with partitions 
        # if type(self.genetic_network)==list:
        #     # if many gene networks - use total_substep_stochastic
            # while delta_t < dt:
            #     # print(f'Elapsed time: {delta_t}')
            #     delta_t += self.total_substep_stochastic(t, dt)
        # else:
        #     for gp in self.gene_products:
        #         if gp.diffusion_rate != 0:
        #             while delta_t < dt:
        #                 #print(f'Elapsed time: {delta_t}')
        #                 delta_t += self.total_substep_stochastic(t, dt)
        #             break
        #     # TODO: do I even need this code?
        #     if delta_t == 0:
        #         self.genetic_network.step_stochastic(t, dt, self.growth_rate)
        '''

    def step(self, t, dt, stochastic=False):
        if self.gene_products and self.biomass:
            for supp,conc in self.supplements.items():
                supp.concentration = conc
            if stochastic:
                self.step_stochastic(t, dt)
            else:
                for s in self.strain:
                    s.genetic_network.step(s.biomass(t), s.growth_rate(t), t, dt)
                # update the exctracellular concentration
                self.external_step(dt)
                self.update_ext_conc(stochastic)



