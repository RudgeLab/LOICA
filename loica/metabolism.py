import numpy as np
from scipy.interpolate import interp1d

def ramp_growth_rate(t, start, slope):
    gr = np.maximum(0, slope*(t-start))
    return(gr)

def ramp_biomass(t, od0, start, slope):
    logod = np.maximum(0, ((t-start)**2)/2)
    od = od0 * np.exp(logod)
    return(od)

def step_growth_rate(t, start):
    gr = 1 * ((t-start)>0)
    return(gr)

def step_biomass(t, od0, start):
    logod = np.maximum(0, t-start)
    od = od0 * np.exp(logod)
    return(od)

def gompertz_growth_rate(t, y0, ymax, um, l):
    A = np.log(ymax/y0)
    gr = um *np.exp((np.exp(1)* um *(l - t))/A - \
            np.exp((np.exp(1)* um *(l - t))/A + 1) + 2)
    return(gr)

def gompertz(t, y0, ymax, um, l):
    A = np.log(ymax/y0)
    log_rel_od = (A*np.exp(-np.exp((((um*np.exp(1))/A)*(l-t))+1)))
    od = y0 * np.exp(log_rel_od)
    return(od)

class Metabolism:
    """
    Context for gene expression, incorporates biomass and growth rate.
    ...
    
    Attributes
    ----------
    name : str, optional
        Name of the metabolism or correponding strain
    """
    def __init__(self, name=None):
        self.name = name

class SimulatedMetabolism(Metabolism):
    """
    Simulated context for gene expression, incorporates biomass and growth rate.
    ...

    Attributes
    ----------
    name : str, optional
        Name of the metabolism or correponding strain
    biomass
        A function of time that describes biomass f(t)=biomass
    growth_rate
        A function of time that describes the growth rate f(t)=growth rate
    """
    def __init__(self, name, biomass, growth_rate):
        super().__init__(name)
        self.biomass = biomass
        self.growth_rate = growth_rate

class DataMetabolism(Metabolism):
    """
    Characterized context for gene expression, incorporates biomass and growth rate.
    ...

    Attributes
    ----------
    name : str, optional
        Name of the metabolism or correponding strain
    fj : Flapjack
        Flapjack instance used to fetch data from
    media : str
        Name of the media to query
    strain : str
        Name of the strain to query
    vector : str
        Name of the vector to query
    biomass_signal : str
        Name of signal to query and use as biomass
    
    Methods
    -------
    biomass(t)
        Return biomass at a given time from characterization data
    growth:rate(t)
        Return growth rate at a given time from characterization data
    """
    def __init__(self, name, fj, media, strain, vector, biomass_signal):
        super().__init__(name)
        gr = fj.analysis(media=media.id, 
                    strain=strain.id,
                    vector=vector.id,
                    type='Expression Rate (indirect)',
                    signal=biomass_signal.id,
                    biomass_signal=biomass_signal.id,
                    pre_smoothing=5,
                    post_smoothing=5).sort_values('Time')
        t = gr.Time.unique()
        gr_prof = gr.groupby('Time').mean().Rate.values
        self.growth_rate_prof = interp1d(t, gr_prof, bounds_error=False, fill_value='extrapolate')

        od = fj.analysis(media=media.id, 
                    strain=strain.id,
                    vector=vector.id,
                    type='Background Correct',
                    signal=biomass_signal.id,
                    biomass_signal=biomass_signal.id
                      ).sort_values('Time')
        t = od.Time.unique()
        biomass_prof = od.groupby('Time').mean().Measurement.values
        self.biomass_prof = interp1d(t, biomass_prof, bounds_error=False, fill_value='extrapolate')

    def biomass(self, t):
        return float(self.biomass_prof(t))

    def growth_rate(self, t):
        return float(self.growth_rate_prof(t))
