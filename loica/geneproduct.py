from random import sample
from .metabolism import convert_to_cells
import numpy as np

class GeneProduct:
    """
    A class that represents a gene product, protein or RNA.

    ...
    
    Attributes
    ----------
    name : str
        Name of the gene product
    init_concentration : int | float
        Initial concentration of the gene product in Molar
    degradation_rate : int | float
        Degradation rate of the gene product
    type_ : str, optional
        Molecular type of the gene product, could be 'PRO' or 'RNA'
    uri : str, optional
        SynBioHub URI
    sbol_comp : SBOL Component, optional
        SBOL Component
    """
    shape = '^'
    def __init__(self, name, init_concentration=0, degradation_rate=0, diffusion_rate=0, uri=None, sbol_comp=None, type_='PRO', color='silver'):
        self.init_concentration = init_concentration
        self.concentration = init_concentration
        self.degradation_rate = degradation_rate
        self.degradation_r = degradation_rate # constant used to reset degradation_rate
        self.diffusion_rate=diffusion_rate
        self.name = name
        self.expression_rate = 0
        self.uri = uri
        self.sbol_comp = sbol_comp
        self.type_ = type_ 
        self.color = color

        self.init_ext_conc = 0
        self.ext_conc = self.init_ext_conc
        self.ext_degr_rate = 0
        self.ext_difference = 0
        self.int_change=0
        self.ext_degraded=0
        self.strain = None

        # test
        self.test = 0


    def initialize(self):
        self.concentration = self.init_concentration
        self.ext_conc = self.init_ext_conc

    def express(self, rate):
        self.expression_rate += rate

    def degrade(self, rate):
        self.degradation_rate += rate

    def step(self, growth_rate, dt, extracellular_volume):
        # this is how much diffused in/out of the cell
        dext_conc_dt = self.diffusion_rate*(self.concentration-self.ext_conc)
        # this section deals with concentration difference between extracellular space and cell 
        # (which is due to different volume)
        if dext_conc_dt<0:
            # convert incoming concentration to moles, then to concentration within cell
            in_moles = dext_conc_dt * extracellular_volume
            converted_conc_change = in_moles / self.strain.cell_volume
            diffusion_cell = converted_conc_change
            diffusion_sample = dext_conc_dt
        elif dext_conc_dt>0:
            # convert outcoming concentration to moles, then to concentration within sample
            in_moles = dext_conc_dt * self.strain.cell_volume
            converted_conc_change = in_moles / extracellular_volume
            diffusion_cell = dext_conc_dt
            diffusion_sample= converted_conc_change
        else:
            diffusion_cell = dext_conc_dt
            diffusion_sample = dext_conc_dt

        # change of concentration within cell
        dconcdt = self.expression_rate - (self.degradation_rate + growth_rate) * self.concentration - diffusion_cell
        
        # test
        without_def = self.expression_rate - (self.degradation_rate + growth_rate) * self.concentration
        self.test = without_def * dt
        
        self.next_concentration = self.concentration + dconcdt * dt
        if self.next_concentration < 0:
            if dext_conc_dt>0:
                # if there is diffusion out of cell, then diffusion and degradation are
                # proportional. This corrects diffusion
                degr_and_gr = ((self.degradation_rate + growth_rate)*self.concentration) * dt
                diff_out = diffusion_cell * dt
                total_minus = degr_and_gr + diff_out
                total_before_minus = self.concentration + self.expression_rate * dt
                proportion_diff = diff_out/total_minus
                corrected_diff = total_before_minus * proportion_diff
                in_moles = corrected_diff * self.strain.cell_volume
                converted_conc_change = in_moles / extracellular_volume
                diffusion_sample= converted_conc_change
            self.concentration = 0
        else:
            self.concentration = self.next_concentration
        # then external concentration change based on number of cells that produce
        # or intake this molecule
        self.ext_difference = diffusion_sample * dt * self.strain.cell_number
        # reset rates
        self.expression_rate = 0
        self.degradation_rate = self.degradation_r
        # test
        # print(f'{self.name} new concentration = {self.concentration}')
        # print(f'{self.name} ext change = {self.ext_difference}')

    def __str__(self):
        return self.name

class Regulator(GeneProduct):
    """
    Representation of a regulatory gene product.
    Child of GeneProduct.
    """
    def __init__(self, name, init_concentration=0, degradation_rate=0, diffusion_rate=0, sbol_comp=None, color='lightgreen'):
        super().__init__(name, init_concentration, degradation_rate, diffusion_rate, sbol_comp,color='lightgreen')
        self.sbol_comp = sbol_comp

class Reporter(GeneProduct):
    """
    Representation of a regulatory gene product.

    Parameters
    ----------
    signal_id : str, optional
        Flapjack ID of the signal that the reporter is associated with.
    color : str, optional
        Color of the reporter
    """
    def __init__(self, name, init_concentration=0, degradation_rate=0, diffusion_rate=0, signal_id=None, color='w', sbol_comp=None):
        super().__init__(name, init_concentration, degradation_rate, diffusion_rate, sbol_comp)
        self.signal_id = signal_id
        self.color = color
        self.sbol_comp = sbol_comp

#Producer
