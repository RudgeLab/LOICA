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
    calibrate(ppod)
        sets particles per OD600, used to calculate cell number
    set_extracel_degr(chemical_name, ext_degr_rate)
        set extracellular degradation rate for GeneProducts with the same name
    set_ext_conc(chemical_name, ext_concentration)
        set starting extracellular concentration for GeneProducts with the same name
    """
    def __init__(self, strain=None, media=None, volume=0.0002):

        self.strain = []
        self.media = media
        self.reporters = []
        self.gene_products = []
        self.growth_rate = []
        self.biomass = []
        self.supplements = {}
        self.volume = volume    # default volume is in L
        self.ppod = 2.66*10**12  # and default ppod is cells/L
        ''' Default ppod(cells per 1 OD600 per volume) has been taken from:
                Yap, P. Y., Trau, D. (2019). 
                DIRECT E.COLI CELL COUNT AT OD600. 
                https://tipbiosystems.com/wp-content/uploads/2020/05/AN102-E.coli-Cell-Count_2019_04_25.pdf
        '''
        self.extracel_vol = volume        

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
                # assign strain to each gene_product
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

    def initialize(self):
        for s in self.strain:
            s.genetic_network.initialize()

    def set_supplement(self, supplement, concentration):
        # TODO: also set how often it is added, how much volume
        self.supplements[supplement] = concentration

    def supplement_is_gp(self, supplement):
        ''' 
            Sets external concentration of gp to supplement concentration is they 
            have the same name
        '''
        # TODO: recalculate concentration based on moles per new volume
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
        ''' 
            Calculates extracellular volume of the sample as well as 
            updates cell_number of each strain
        '''
        extracel_v = self.extracel_vol
        for s in self.strain:
            current_cell_n = convert_to_cells(s.biomass(t), self.ppod, self.volume)
            difference = s.cell_number - current_cell_n 
            s.cell_number = current_cell_n
            # if there are more cells, difference is negative, extracellular volume 
            # decreases
            extracel_v += difference * s.cell_volume
        # update external concentration due to volume change:
        if t != 0:
            for group in self.gene_products:
                moles = group[0].ext_conc / self.extracel_vol
                updated_ext_conc = moles / extracel_v
                for gp in group: 
                    gp.ext_conc = updated_ext_conc
        self.extracel_vol = extracel_v

    
    def external_step(self, dt):
        """ 
            Calculates the change in the extracellular concentration
            due to degradation
        """
        for group in self.gene_products: 
            if group[0].ext_degr_rate != 0:
                ext_degr = group[0].ext_conc * group[0].ext_degr_rate * dt
                for gp in group:
                    gp.ext_degraded = ext_degr

    def update_ext_conc(self, t):  
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
                new_ext_conc = self.catch_negative_conc(group)
            for gp in group:
                gp.ext_conc = new_ext_conc

    def catch_negative_conc(self, group):
        ''' 
            used if external concentration becomes negative due to multiple strains/
            multiple bacteria taking molecules out of the extracellular space 
            simultaneously.
            
            Diffusion is recalculated to be proportionate and molecules are "returned" 
            from the cell.
        '''
        new_ext_conc = group[0].ext_conc
        # list of gene products that diffuse out of the extracellular space
        diffused_out = []
        for gp in group:
            if gp.ext_difference<0:
                diffused_out.append(gp)
            else:
                new_ext_conc += gp.ext_difference
        # calculate what would be the change of concentration ideally
        ideal_minus = 0
        for gp in diffused_out:
            ideal_minus -= gp.ext_difference
        if group[0].ext_degraded > 0:
            ideal_minus += group[0].ext_degraded
        for gp in diffused_out:
            # calculate proportion of concentration each strain would take
            # then multiply by the available concentration to get the concentration
            # change that happened
            can_take = ((-gp.ext_difference)/ideal_minus*new_ext_conc)
            # calculate "extra" concentration each cell has taken
            extra_taken = (-gp.ext_difference-can_take)/gp.strain.cell_number 
            # convert concentration to concentration within cell:
            in_moles = extra_taken * self.extracel_vol
            extra_conc_converted = in_moles / gp.strain.cell_volume
            # correct the internal gp concentration
            fixed_conc = gp.concentration - extra_conc_converted
            if fixed_conc < 0:
                # this might happen if extra concentration that diffused into the cell was 
                # degraded straight away
                gp.concentration = 0
            else:
                gp.concentration = fixed_conc
        return 0
                
    def step(self, t, dt):
        ''' deterministic '''
        if self.gene_products and self.biomass:
            for supp,conc in self.supplements.items():
                # TODO: change
                supp.concentration = conc
                self.supplement_is_gp(supp)
            # calculate cell number, extracellular volume and external concentration
            self.extracel_volume(t)
            # step
            for s in self.strain:
                s.genetic_network.step(s.growth_rate(t), t, dt, self.extracel_vol)
            # update the exctracellular concentration
            self.external_step(dt)
            self.update_ext_conc(t)




