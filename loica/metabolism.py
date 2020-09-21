import numpy as np

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
        self.biomass_func = biomass_func
        self.growth_rate_func = growth_rate_func

    def biomass(self, t):
        return self.biomass_func(t)

    def growth_rate(self, t):
        return self.growth_rate_func(t)

