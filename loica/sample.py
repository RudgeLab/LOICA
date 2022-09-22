from itertools import groupby
from .metabolism import convert_to_cells
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
        strains which are in the sample
    media : str
        Name of the media in the sample
    volume: int | float
        Sample volume. By default is set to 1.36E-08 L to represent an 85*100*1600 um trap
        in microfluidic device from https://doi.org/10.1126%2Fscience.aaa3794
    
     Methods
    -------
    add_supplement(supplement, concentration)
        establishes the concentration of Supplement
    set_extracel_degr(chemical_name, ext_degr_rate)
        set extracellular degradation rate for GeneProducts with the same name
    set_ext_conc(chemical_name, ext_concentration)
        set starting extracellular concentration for GeneProducts with the same name
    """
    def __init__(self, strain=None, media=None, volume=1.36*10**-8):

        self.strain = []
        self.media = media
        self.reporters = []
        self.gene_products = []
        self.growth_rate = []
        self.biomass = []
        self.supplements = {}
        self.volume = volume    # default volume is in L
        self.ppod = 3.26*10**11 # and default ppod is cells/L
        ''' Default ppod(cells per 1 OD600 per volume) has been taken from example
            data for iGEM 2019 calibration protocol:
            example data - https://static.igem.org/mediawiki/2019/2/29/IGEM_2019_Data_Analysis_Template_Particle_Standard_Curve_Protocol_v2.xlsx
            protocol - Beal, J., Haddock-Angelli, T., Gershater, M., Sanchania, V., Buckley-Taylor, R., Baldwin, G., Farny, N., Tennant, R., &#38; Rutten, P. (2020). 
                       Calibration Protocol - Plate Reader Abs600 (OD) Calibration with Microsphere Particles V.3. https://www.protocols.io/view/calibration-protocol-plate-reader-abs600-od-calibr-ewov1ypoovr2/v3
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
        ''' sets particle per OD600 - used to convert absorbance to cell number
            assumes that all strains have same ppod value'''
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
        if t!=0 and extracel_v!=self.extracel_vol:
            for group in self.gene_products:
                moles = group[0].ext_conc * self.extracel_vol
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
                raise ValueError(f'Negative {group[0].name} extracellular concentration')
            for gp in group:
                gp.ext_conc = new_ext_conc
                
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





