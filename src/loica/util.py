import pickle
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import least_squares

def forward_model_growth(
    Dt=0.05,
    sim_steps=10,
    muval=[0]*100,
    od0=0,
    nt=100
):
    od_list, t_list = [],[]
    od = od0
    for t in range(nt):
        od_list.append(od)
        t_list.append([t * Dt])
        mu = muval[t]
        for tt in range(sim_steps):
            doddt = mu * od
            nextod = od + doddt * Dt/sim_steps
            od = nextod


    aod = np.array(od_list).transpose()
    tt = np.array(t_list).transpose()
    return aod,tt


def residuals_growth(data, epsilon, dt, t, n_gaussians): 
    def func(x): 
        od0 = x[0]
        muval = np.zeros_like(t)
        means = np.linspace(t.min(), t.max(), n_gaussians)
        vars = [(t.max()-t.min())/n_gaussians] * n_gaussians 
        heights = x[1:]
        for mean,var,height in zip(means, vars, heights):
            gaussian = height * np.exp(-(t-mean)*(t-mean) / var / 2) / np.sqrt(2 * np.pi * var)
            muval = muval + gaussian

        od,tt = forward_model_growth(
                    Dt=dt,
                    muval=muval,
                    od0=od0,
                    nt=len(t)
                )
        model = od
        residual = (data - model)  # / tt.ravel()[1:]
        tikhonov = heights
        result = np.concatenate((residual, epsilon * tikhonov))
        return result
    return func


def characterize_growth(flapjack, 
        vector, 
        media, 
        strain, 
        biomass_signal, 
        n_gaussians, 
        epsilon
        ):
    # Characterize growth rate profile
    biomass = flapjack.analysis(type='Background Correct', 
                        vector=vector,
                        media=media,
                        strain=strain,
                        signal=biomass_signal,
                        biomass_signal=biomass_signal
                     )
    biomass = biomass.sort_values('Time')
    t = biomass.Time.unique()
    dt = np.mean(np.diff(t))
    nt = len(t)

    biomass = biomass.groupby('Time').mean().Measurement.values
    lower_bounds = [0] + [0]*n_gaussians
    upper_bounds = [100] + [50]*n_gaussians
    bounds = [lower_bounds, upper_bounds]

    data = biomass
    res = least_squares(
            residuals_growth(data, epsilon=epsilon, dt=dt, t=t, n_gaussians=n_gaussians), 
            [0.01] + [1]*n_gaussians, 
            bounds=bounds
            )
    init_biomass = res.x[0]
    profile = np.zeros_like(t)
    means = np.linspace(t.min(), t.max(), n_gaussians)
    vars = [(t.max()-t.min())/n_gaussians] * n_gaussians 
    heights = res.x[1:]
    for mean,var,height in zip(means, vars, heights):
        gaussian = height * np.exp(-(t-mean)*(t-mean) / var / 2) / np.sqrt(2 * np.pi * var)
        profile = profile + gaussian
    mu_profile = interp1d(t, profile, fill_value='extrapolate', bounds_error=False)

    od,tt = forward_model_growth(
                    Dt=dt,
                    muval=profile,
                    od0=init_biomass,
                    nt=nt
                )
    return init_biomass, od, mu_profile

def save_loica(obj, filename):
    with open(f'{filename}.obj', 'wb') as outp:
        pickle.dump(obj, outp, pickle.HIGHEST_PROTOCOL)

def load_loica(filename):
    with open(f'{filename}.obj', 'rb') as inp:
        loica = pickle.load(inp)
    return loica