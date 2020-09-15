import pandas as pd
import numpy as np
from . import genetic_network, metabolism



class Assay:
    def __init__(self, samples, n_measurements, interval):
        '''
        Assay measures a set of samples in parallel at a set of timepoints
        Connects to flapjack to generate data, and to fit parameters to data

        '''
        self.samples = samples
        self.n_measurements = n_measurements
        self.interval = interval
        self.measurements = pd.DataFrame()

    def run(self, substeps=10):
        '''
        Run the assay measuring at specified time points, with simulation time step dt
        '''
        #substeps = self.interval / dt
        dt = self.interval / substeps
        for sample_id, sample in enumerate(self.samples):
            # Integrate models
            for t in range(self.n_measurements):
                for tt in range(substeps):
                    time = t * self.interval + tt * dt
                    sample.step(time, dt)
                # Record measurements of fluorescence
                signals = sample.signals
                for name, sig in signals.items():
                    row = {'Time': time, 'Signal': name, 'Measurement': sig * sample.biomass, 'Sample':sample_id}
                    self.measurements = self.measurements.append(row, ignore_index=True)
                # Record measurement of biomass
                row = {'Time': time, 'Measurement': sample.biomass, 'Signal':'Biomass', 'Sample':sample_id}
                self.measurements = self.measurements.append(row, ignore_index=True)
                
    def upload(self, 
            flapjack, 
            name,
            temp,
            description,
            study,
            media,
            strain,
            signals
            ):
        assay = flapjack.create('assay', name=name, study=study, temperature=temp, description=description)
        for sample_id, sample in enumerate(samples):
            sample = fj.create('sample',
                                row=sample_id, col=1,
                                media=media.id[0],
                                strain=strain.id[0],
                                vector=vector.id[0],
                                assay=assay.id[0],
                                )

