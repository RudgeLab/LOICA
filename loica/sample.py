from itertools import groupby
from .metabolism import convert_to_cells
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
    volume: int | float
        Sample volume. By default is set to 0.0002 L (200 ul) to represent a well in 96-well plate
    
     Methods
    -------
    add_supplement(supplement, concentration)
        establishes the concentration of Supplement
    """
    def __init__(self, strain=None, media=None, volume=0.0002): # resources=None

        self.strain = []
        self.media = media
        self.reporters = []
        self.gene_products = []
        self.growth_rate = []
        self.biomass = []
        self.supplements = {}
        self.volume = volume    # default volume is in L
        self.ppod = 2.66*10**6  # and default ppod is cells/L
        ''' Default ppod(cells per 1 OD600 per volume) has been taken from:
                Yap, P. Y., Trau, D. (2019). 
                DIRECT E.COLI CELL COUNT AT OD600. 
                https://tipbiosystems.com/wp-content/uploads/2020/05/AN102-E.coli-Cell-Count_2019_04_25.pdf
        '''
        self.extracel_vol = volume
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

        # create an options list for semi-stochastic simulation 
        # (used in self.total_substep_stochastic(self))
        self.options = self.strain[:]
    
    def calibrate(self, ppod):
        ''' sets particle per OD600 - used to convert absorbance to cell number'''
        self.ppod = ppod

    def set_extracel_degr(self, chemical_name, ext_degr_rate):
        ''' 
            this method ensures that all gene products with the same identity 
            have the same extracellular degradation rate
        '''
        for group in self.gene_products:
            if group[0].name == chemical_name:
                for gp in group:
                    gp.ext_degr_rate = ext_degr_rate
        
        # this ensures that self.st_external_substep() is called later in self.total_substep_stochastic()
        if "extracellular space" not in self.options:
            self.options.append("extracellular space")

    def initialize(self):
        for s in self.strain:
            s.genetic_network.initialize()

    def set_supplement(self, supplement, concentration):
        self.supplements[supplement] = concentration

    def supplement_is_gp(self, supplement):
        ''' set external concentration of gp to supplement concentration is they have the same name'''
        for group in self.gene_products:
            if group[0].name == supplement.name:
                for gp in group:
                    gp.ext_conc += supplement.concentration

    def set_ext_conc(self, chemical_name, ext_concentration):
        ''' set initial external concentration '''
        
        for group in self.gene_products:
            if group[0].name == chemical_name:
                for gp in group:
                    gp.init_ext_conc = ext_concentration

    def set_regulator(self, name, strain, concentration):
        for s in self.strain:
            if strain == s.name:
                for reg in s.regulators:
                    if reg.name == name:
                        reg.init_concentration = concentration
                    else: pass

    def set_reporter(self, name, strain, concentration):
        for s in self.strain:
            if strain == s.name:
                for rep in s.regulators:
                    if rep.name == name:
                        rep.init_concentration = concentration
                    else: pass

    def extracel_volume(self, t):
        ''' this function calculates extracellular volume of the sample '''
        extracel_v = self.extracel_vol
        for s in self.strain:
            s.cell_number = convert_to_cells(s.biomass(t), self.ppod, self.volume)
            extracel_v -= s.cell_number * s.cell_volume
        self.extracel_vol = extracel_v

    
    def external_step(self, dt):
        """ 
            method that calculates the change in the extracellular concentration
            due to degradation

            deterministic
        """
        for group in self.gene_products: 
            ext_degr = group[0].ext_conc * group[0].ext_degr_rate * dt
            for gp in group:
                gp.ext_degraded = ext_degr

    def update_ext_conc(self, t, stochastic=False):  
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
                new_ext_conc = self.catch_negative_conc(group, t, stochastic)
            for gp in group:
                gp.ext_conc = new_ext_conc
            # test 
            if t<5:
                for gp in group:
                    print(f'''After update:
                    {gp.name} new ext conc = {gp.ext_conc}
                    int conc in {gp.strain.name} = {gp.concentration}''')

    def catch_negative_conc(self, group, t, stochastic=False):
        ''' 
            used if external concentration becomes negative due to multiple strains
            taking molecules out of the extracellular space simultaneously.
            randomly updates extracellular concentration
            if result is negative, diffusion does not happen and molecules
            are returned to the cell
        '''
        new_ext_conc = group[0].ext_conc
        if stochastic:
            # only needed if using for semi-stochastic and fully stochastic with partitions 
            # needs to be revised
            ext_options = ["degradation"]
            for gp in group:
                ext_options.append(gp) 
            shuffle(ext_options)
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
            # test
            if t<5:
                print("triggered negative update")
            # list of gene products that diffuse out of the extracellular space
            diffused_out = []
            for gp in group:
                if gp.ext_difference<0:
                    diffused_out.append(gp)
                else:
                    new_ext_conc += gp.ext_difference
            # how much would ideally be taken out and degraded
            ideal_minus = 0
            for gp in diffused_out:
                # calculate what would be the change of concentration if there were 
                # enough resources
                ideal_minus -= gp.ext_difference
            if group[0].ext_degraded > 0:
                ideal_minus += group[0].ext_degraded
            # test
            if t<5:
                print(f''' {group[0].name} ideal diffusion out with degr = {ideal_minus}
                while ext conc after addition to it = {new_ext_conc}''')
            for gp in diffused_out:
                # calculate proportion of molecules each strain would take in from total
                # then multiply by the available concentration to get the concentration
                # change that happened
                can_take = ((-gp.ext_difference)/ideal_minus*new_ext_conc)
                # calculate "extra" concentration of the gp in the strain it belongs to
                extra_taken = (-gp.ext_difference-can_take)/gp.strain.cell_number 
                # convert concentration to concentration within cell:
                in_moles = extra_taken * self.extracel_vol
                extra_conc_converted = in_moles / gp.strain.cell_volume
                # correct the internal gp concentration
                fixed_conc = gp.concentration - extra_conc_converted
                # test
                if t<5:
                    print(f''' {gp.strain.name} can take = {can_take}
                    cell_number = {gp.strain.cell_number}
                    it has taken in extra {extra_taken} per cell
                    which is {extra_conc_converted} within cell
                    conc before correction = {gp.concentration}
                    after correction = {fixed_conc}
                    ''')
                if fixed_conc < 0:
                    # this might happen if extra concentration that diffused into the cell was 
                    # degraded straight away
                    gp.concentration = 0
                else:
                    gp.concentration = fixed_conc
                
            return 0
                
    def st_external_substep(self, tau_=None):
        """ 
            method similar to GeneticNetwork.substep()
            but only for extracellular degradation. Used in simulation with partitions.

            stochastic
        """
        # Propensities
        a = []

        # Compute propensities for degradation of gene products
        for group in self.gene_products:
            # degradation reaction
            a.append(group[0].ext_degr_rate * group[0].ext_conc)
            # reset how much degraded
            group[0].ext_degraded = 0
        
        # Make list of propensities into array
        a = np.array(a)
        # Total of propensities
        A = a.sum()

        if A == 0:
            if tau_:
                tau = tau_
            else:
                tau = 0
            return tau
        
        # Time step
        if tau_ and tau_ != 0:
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

        # Return elapsed time
        return tau
        
    def total_substep_stochastic(self, t=0, dt=0.1, type='stochastic'):
        '''
            method that links stochastic substeps for each genetic network and 
            extracellular space.
            semi-stochastic or fully stochastic
            with partitions
        '''
        # shuffle the list
        shuffle(self.options)

        # get tau by running substep for the first item in the shuffled list
        # and then use this tau in other substeps
        if self.options[0]=="extracellular space":
            tau = self.st_external_substep()
            for s in self.options[1:]:
                if type=='full+comp':               
                    new_tau = s.genetic_network.substep_stochastic(t, dt, s.growth_rate(t), s.biomass(t), tau)
                elif type=='semi+comp':
                    new_tau = s.genetic_network.substep_semistochastic(t, dt, s.growth_rate(t), s.biomass(t), tau)
                    tau = new_tau
            self.update_ext_conc(stochastic=True)
        else:
            if type=='full+comp':
                tau = self.options[0].genetic_network.substep_stochastic(t, dt, self.options[0].growth_rate(t), self.options[0].biomass(t))
            elif type=='semi+comp':
                tau = self.options[0].genetic_network.substep_semistochastic(t, dt, self.options[0].growth_rate(t), self.options[0].biomass(t))
            for o in self.options[1:]:
                if o == "extracellular space":
                    new_tau = self.st_external_substep(tau_=tau)
                    tau = new_tau
                else:
                    if type=='full+comp':
                        new_tau = o.genetic_network.substep_stochastic(t, dt, o.growth_rate(t), o.biomass(t), tau)
                    elif type=='semi+comp':
                        new_tau = o.genetic_network.substep_semistochastic(t, dt, o.growth_rate(t), o.biomass(t), tau)
                        tau = new_tau
            self.update_ext_conc(stochastic=True)
                    
        return tau

    def substep_stochastic(self, t=0, dt=0.1, ppod=2.66*10**9):
        ''' fully stochasic, only one reaction happens out of all'''
        # Propensities
        a = []

        # Compute expression rates
        for s in self.strain:
            for op in s.genetic_network.operators:
                expression_rate = op.expression_rate(t, dt)
                output = op.output
                if type(op.output)!=list:
                    output = [output]
                for o in output:
                    o.express(expression_rate)

        # Compute propensities for production, degradation, diffusion in and out of gene products
        for group in self.gene_products:
            for gp in group:
                # Production reaction
                a.append(gp.expression_rate)
                # Degradation
                a.append((gp.degradation_rate + gp.strain.growth_rate(t)) * gp.concentration)
                # Diffusion out of cell
                a.append(gp.diffusion_rate*gp.concentration)
                # Difusion into the cell
                a.append(gp.diffusion_rate*gp.ext_conc)
                # External degradation
                a.append(gp.ext_conc*gp.ext_degr_rate)
                # reset values
                gp.ext_difference = 0

        # Make list of propensities into array
        a = np.array(a)
        # Total of propensities
        A = a.sum()
        
        # Time step
        tau = 1/A * np.log(1/np.random.random())
        # Random number to select next reaction
        a_i = np.random.random() * A

        # Find reaction and update gene product levels
        x = 0
        complete = False
        for group in self.gene_products:
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
                    gp.ext_difference = convert_to_cells(gp.strain.biomass(t), ppod, self.volume)
                    complete = True
                    break
                elif a_i < np.sum(a[:i*5+4]):
                    # Diffusion of geneproduct gp into the cell
                    gp.concentration += 1
                    gp.ext_difference = - convert_to_cells(gp.strain.biomass(t), ppod, self.volume)
                    complete = True
                    break
                elif a_i < np.sum(a[:i*5+5]):
                    # External degradation of geneproduct gp
                    gp.ext_difference -= 1
                    complete = True
                    break
            if complete:
                break
            x += len(group)
        
        # Reset expression rates for next step and update external concentration
        # TODO: catch negative external concentration
        for group in self.gene_products:
            concentration_change = 0
            for gp in group:
                gp.expression_rate = 0
                concentration_change += gp.ext_difference   
            new_ext_conc = group[0].ext_conc + concentration_change
            for gp in group:
                gp.ext_conc = new_ext_conc
        
        # Return elapsed time
        return tau

    def step_stochastic(self, t=0, dt=0.1, type='full_stochastic', ppod=2.66*10**9):
        ''' 
            similar to GeneticNetwork.step_stochastic(). 
            Use but uses either 
            self.total_substep_stochastic() or GeneticNetwork.step_stochastic()
        '''
        delta_t = 0
        if type=='full_stochastic':
            while delta_t < dt:
                # print(f'Elapsed time: {delta_t}')
                delta_t += self.substep_stochastic(t, dt, ppod=ppod)
        elif type=='semi+comp' or type=='full+comp':
            while delta_t < dt:
                # print(f'Elapsed time: {delta_t}')
                delta_t += self.total_substep_stochastic(t, dt, type)
                if delta_t == 0:
                    return

    def step(self, t, dt, stochastic=False):
        if self.gene_products and self.biomass:
            for supp,conc in self.supplements.items():
                # TODO: change
                supp.concentration = conc
                self.supplement_is_gp(supp)
            if stochastic:
                if type(stochastic)==str:
                    self.step_stochastic(t, dt, type=stochastic)
                else:
                    self.step_stochastic(t, dt, self.ppod)
            else:
                self.extracel_volume(t)
                for s in self.strain:
                    s.genetic_network.step(s.growth_rate(t), t, dt, self.extracel_vol)
                # update the exctracellular concentration
                self.external_step(dt)
                # test
                if t<5:
                    for group in self.gene_products:
                        for gp in group:
                            if gp.ext_degr_rate > 0:
                                print(f'''After degradation 
                                {gp.name} in {gp.strain.name} ext conc = {gp.ext_conc}''')
                self.update_ext_conc(t, stochastic)



