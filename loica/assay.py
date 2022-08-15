import pandas as pd
import numpy as np
from tqdm import tqdm
from . import genetic_network, metabolism



class Assay:
    """
    Assay measures a set of samples in parallel at a set of timepoints.
    Connects to flapjack to generate data, and to fit parameters to data.
    
    ...
    
    Attributes
    ----------
    samples : List[Sample]
        List of Samples that belongs to the Assay
    n_measurements : int
        Number of measurements to take
    interval : int
        Time in hours between each measurements
    name : str
        Name of the Assay
    description: str
        Descriptioin of the Assay
    biomass_signal_id : int
        Flapjack ID of the Assay that is associated with the Assay

    Methods
    -------
    run(substeps=10, nsr=0, biomass_bg=0, fluo_bg=0, mode=None, ppod=None)
        Runs the Assay time series
        mode can be either:
        (empty) - gives overall fluorescence and biomass measurements
        'track all' - raw numbers from the simulation. Records all GeneProducts in each
                      strain and extracellular space. Also records biomass of each strain.
        'single_cell' - fluorescence measurements of a cell per strain
    upload(flapjack, study)
        Upload the data produced by running the Assay to Flapjack into the Study
    """
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


    def run(self, substeps=10, nsr=0, biomass_bg=0, fluo_bg=0, mode=None, ppod=None):
        '''
        Run the assay measuring at specified time points, with simulation time step dt
        '''
        #substeps = self.interval / dt
        dt = self.interval / substeps
        n_samples = len(self.samples)
        with tqdm(total=100) as pbar:
            for sample_id, sample in enumerate(self.samples):
                sample.initialize()
                # Calibrate sample so cell number could be calculated from absorbance
                if ppod:
                    sample.calibrate(ppod)
                # Integrate models
                for t in range(self.n_measurements):
                    # print(f'Current t={t}')
                    time = t * self.interval
                    if mode == 'track_all':
                    # Record raw measurements of each gene product both intracellularly 
                    # and extracellularly
                        for tt in range(substeps):
                            for group in sample.gene_products:
                                row = {
                                        'Time': time, 
                                        'Signal_id': None, 
                                        'Signal': f'Extracellular {group[0].name}', 
                                        'Measurement': group[0].ext_conc,
                                        'Sample':sample_id
                                        }
                                self.measurements = self.measurements.append(row, ignore_index=True)
                                for gp in group:
                                    row = {
                                        'Time': time, 
                                        'Signal_id': None, 
                                        'Signal': f'{gp.name} in {gp.strain.name}', 
                                        'Measurement': gp.concentration,
                                        'Sample':sample_id
                                        }
                                    self.measurements = self.measurements.append(row, ignore_index=True)
                                   
                            total = 0
                            # biomass of each strain recorded separately
                            for strain in sample.strain:
                                row = {
                                        'Time': time, 
                                        'Signal_id': None, 
                                        'Measurement': strain.biomass(time),
                                        'Signal':f'{strain.name} biomass', 
                                        'Sample':sample_id
                                        }
                                self.measurements = self.measurements.append(row, ignore_index=True)
                                total += strain.biomass(time)
                            # total biomass
                            row = {
                                    'Time': time, 
                                    'Signal_id': None, 
                                    'Measurement': total,
                                    'Signal':'Total biomass', 
                                    'Sample':sample_id
                                    }
                            self.measurements = self.measurements.append(row, ignore_index=True)
                            time = t * self.interval + tt * dt
                            sample.step(time, dt)

                    elif mode=='single_cell':
                        for strain in sample.strain:
                            for reporter in strain.reporters:
                                sig = reporter.concentration
                                signal_id = reporter.signal_id
                                signal_name = reporter.name
                                noise = np.random.normal(scale=np.sqrt(nsr))
                                meas = fluo_bg + reporter.ext_conc 
                                noisy_meas = (1 + noise) * meas
                                noise_bg = np.random.normal(scale=np.sqrt(nsr))
                                corr_meas = noisy_meas - (1 + noise_bg) * fluo_bg
                                row = {
                                        'Time': time, 
                                        'Signal_id': signal_id, 
                                        'Signal': f'{signal_name} in {strain.name}', 
                                        'Measurement': corr_meas,
                                        'Sample':sample_id
                                        }
                                self.measurements = self.measurements.append(row, ignore_index=True)
                    else:
                        for reporter in sample.reporters:
                            signal_id = reporter.signal_id
                            signal_name = reporter.name
                            noise = np.random.normal(scale=np.sqrt(nsr))
                            meas = fluo_bg + reporter.ext_conc 
                            for r in sample.reporters:
                                if r.name == reporter.name:
                                    sig = r.concentration
                                    # multiply signal by the number of cells showing it
                                    meas += sig * reporter.strain.cell_number
                            noisy_meas = (1 + noise) * meas
                            noise_bg = np.random.normal(scale=np.sqrt(nsr))
                            corr_meas = noisy_meas - (1 + noise_bg) * fluo_bg
                            row = {
                                    'Time': time, 
                                    'Signal_id': signal_id, 
                                    'Signal':signal_name, 
                                    'Measurement': corr_meas,
                                    'Sample':sample_id
                                    }
                            self.measurements = self.measurements.append(row, ignore_index=True)

                        # Record measurement of biomass
                        noise = np.random.normal(scale=np.sqrt(nsr))
                        noise_bg = np.random.normal(scale=np.sqrt(nsr))
                        meas = biomass_bg
                        for strain in sample.strain:
                            meas += strain.biomass(time)
                        noisy_meas = (1 + noise) * meas
                        corr_meas = noisy_meas - (1 + noise_bg) * biomass_bg
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
                pbar.update(1 / n_samples * 100)
            pbar.close()
                
    def upload(self, flapjack, study):
        assay = flapjack.create('assay', 
                                    name=self.name, 
                                    study=study, 
                                    temperature=0, 
                                    machine='Loica',
                                    description=self.description)
        n_samples = len(self.samples)
        with tqdm(total=100) as pbar:
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
                pbar.update(1 / n_samples * 100)
            pbar.close()
