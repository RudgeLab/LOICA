import numpy as np
from scipy.interpolate import interp1d

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
    def __init__(self):
        pass

class SimulatedMetabolism(Metabolism):
    def __init__(self, biomass_func, growth_rate_func):
        super().__init__()
        self.biomass = biomass_func
        self.growth_rate = growth_rate_func

class DataMetabolism(Metabolism):
    def __init__(self, fj, media, strain, vector, biomass_signal):
        super().__init__()
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
