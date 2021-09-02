import pandas as pd
import numpy as np
from . import genetic_network, metabolism



class Assay:
    def __init__(self, 
            samples, 
            n_measurements, 
            interval,
            name='Loica assay',
            description='',
            biomass_signal_id=None
            ):
        '''
        Assay measures a set of samples in parallel at a set of timepoints
        Connects to flapjack to generate data, and to fit parameters to data

        '''
        self.samples = samples
        self.n_measurements = n_measurements
        self.interval = interval
        self.measurements = pd.DataFrame()
        self.name = name
        self.description = description
        self.biomass_signal_id = biomass_signal_id

    def run(self, substeps=10, noise_std=0, biomass_bg=0, fluo_bg=0):
        '''
        Run the assay measuring at specified time points, with simulation time step dt
        '''
        #substeps = self.interval / dt
        dt = self.interval / substeps
        for sample_id, sample in enumerate(self.samples):
            sample.initialize()
            # Integrate models
            for t in range(self.n_measurements):
                time = t * self.interval
                # Record measurements of fluorescence
                for reporter in sample.reporters:
                    sig = reporter.concentration
                    signal_id = reporter.signal_id
                    signal_name = reporter.name
                    noise_val = np.random.normal(scale=noise_std)
                    meas = sig * sample.biomass(time) + fluo_bg
                    noisy_meas = meas + noise_val
                    corr_meas = noisy_meas - fluo_bg
                    row = {
                            'Time': time, 
                            'Signal_id': signal_id, 
                            'Signal':signal_name, 
                            'Measurement': corr_meas,
                            'Sample':sample_id
                            }
                    self.measurements = self.measurements.append(row, ignore_index=True)
                # Record measurement of biomass
                noise_val = np.random.normal(scale=noise_std)
                meas = sample.biomass(time) + biomass_bg
                noisy_meas = meas + noise_val
                corr_meas = noisy_meas - biomass_bg
                row = {
                        'Time': time, 
                        'Signal_id': self.biomass_signal_id, 
                        'Measurement': corr_meas,
                        'Signal':'Biomass', 
                        'Sample':sample_id
                        }
                self.measurements = self.measurements.append(row, ignore_index=True)
                # Compute next time step
                for tt in range(substeps):
                    time = t * self.interval + tt * dt
                    sample.step(time, dt)
                
    def upload(self, flapjack, study):
        assay = flapjack.create('assay', 
                                    name=self.name, 
                                    study=study, 
                                    temperature=0, 
                                    machine='Loica',
                                    description=self.description)
        for sample_id,sample in enumerate(self.samples):
            fj_sample = flapjack.create('sample',
                                row=sample_id, col=1,
                                media=sample.media,
                                strain=sample.strain,
                                vector=sample.vector,
                                assay=assay.id[0],
                                )

            supplements = []
            for supp,conc in sample.supplements.items():
                chemical = flapjack.get('chemical', name=supp.name)
                if len(chemical)==0:
                    chemical = flapjack.create('chemical', name=supp.name, description='Testing')
                name = supp.name + f' {conc}'
                s = flapjack.get('supplement', name=name, concentration=conc, chemical=chemical.id[0])
                if len(s)==0:
                    s = flapjack.create('supplement', name=name, concentration=conc, chemical=chemical.id[0])
                supplements.append(s.id[0])
            if len(supplements):
                flapjack.patch('sample', fj_sample.id[0], supplements=supplements)

            meas = self.measurements[self.measurements.Sample==sample_id]
            for signal_id,data in meas.groupby('Signal_id'):
                if signal_id:
                    flapjack.upload_measurements(data, signal=[signal_id], sample=fj_sample.id)

