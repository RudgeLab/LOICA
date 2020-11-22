import numpy as np
from numpy.fft import fft, ifft, fftfreq
from scipy.optimize import least_squares
from scipy.interpolate import interp1d

class Source:
    def __init__(self, output, rate, profile=None):
        if profile:
            self.profile = profile
        else:
            def profile(t):
                return 1
            self.profile = profile
        self.rate = rate
        self.output = output

    def expression_rate(self, t, dt):
        return self.rate * self.profile(t)

    def forward_model(
        self,
        Dt=0.25,
        sim_steps=10,
        odval=[1]*97,
        profile=[1]*97,
        gamma=0,
        p0=0,
        nt=100
    ):
        p1_list,od_list, A_list,t_list = [],[],[],[]
        p1 = p0
        for t in range(nt):
            p1_list.append(p1)
            t_list.append([t * Dt])
            od = odval[t]
            tt = t*Dt
            prof = profile[t]
            for tt in range(sim_steps):
                nextp1 = p1 + (odval[t]*profile[t] - gamma*p1) * Dt / sim_steps
                p1 = nextp1


        ap1 = np.array(p1_list).transpose()
        tt = np.array(t_list).transpose()
        t = np.arange(nt) * Dt
        return ap1,tt

    def residuals(self, data, p0, odval, dt, t, n_gaussians, epsilon): 
        def func(x): 
            nt = len(t)
            means = np.linspace(t.min(), t.max(), n_gaussians)
            vars = [(t.max()-t.min())/n_gaussians]*n_gaussians 
            p0 = x[0]
            gamma = x[1]
            heights = x[2:]
            profile = np.zeros_like(t)
            for mean,var,height in zip(means, vars, heights):
                gaussian = height * np.exp(-(t-mean)*(t-mean) / var / 2) / np.sqrt(2 * np.pi * var)
                profile = profile + gaussian
            p,tt = self.forward_model(
                        Dt=dt,
                        odval=odval,
                        profile=profile,
                        nt=nt,
                        p0=p0,
                        gamma=gamma
                    )
            model = p[1:]
            tikhonov = heights * epsilon
            #tikhonov = np.diff(profile) * epsilon
            ntimes = len(t)*dt - tt.ravel()[1:]
            residual = (data[1:] - model)  # / tt.ravel()[1:] 
            return np.concatenate((residual, tikhonov))
        return func

    def characterize(self, flapjack, vector, media, strain, signal, biomass_signal, n_gaussians, epsilon):
        expression_df = flapjack.analysis(media=media, 
                            strain=strain,
                            vector=vector,
                            signal=signal,
                            type='Background Correct',
                            biomass_signal=biomass_signal
                            ).sort_values(['Sample', 'Time'])
        t = expression_df.groupby('Time').mean().index.values
        dt = np.diff(t).mean()
        expression = expression_df.groupby('Time').mean().Measurement.values

        biomass_df = flapjack.analysis(media=media, 
                            strain=strain,
                            vector=vector,
                            signal=biomass_signal,
                            type='Background Correct',
                            biomass_signal=biomass_signal
                            ).sort_values(['Sample', 'Time'])
        biomass = biomass_df.groupby('Time').mean().Measurement.values

        nt = len(t)

        # Bounds for fitting
        lower_bounds = [0] + [0] + [0]*n_gaussians
        upper_bounds = [1e8] + [5] + [1e8]*n_gaussians
        bounds = [lower_bounds, upper_bounds]
        '''
            gamma = x[0]
            profile = x[1:]
        '''
        data = expression.ravel()
        self.residuals_func0 = self.residuals(
                    data, data[0], biomass, epsilon=0, dt=dt, t=t, n_gaussians=n_gaussians
                    )
        self.residuals_func = self.residuals(
                    data, data[0], biomass, epsilon=epsilon, dt=dt, t=t, n_gaussians=n_gaussians
                    )
        res = least_squares(
                self.residuals_func, 
                [0] + [0] + [100]*n_gaussians, 
                bounds=bounds
                )
        self.res = res

        self.p0 = res.x[0]
        self.gamma = res.x[1]

        profile = np.zeros_like(t)
        means = np.linspace(t.min(), t.max(), n_gaussians)
        vars = [(t.max()-t.min())/n_gaussians] * n_gaussians 
        heights = res.x[2:]
        for mean,var,height in zip(means, vars, heights):
            gaussian = height * np.exp(-(t-mean)*(t-mean) / var / 2) / np.sqrt(2 * np.pi * var)
            profile = profile + gaussian
        self.rate = profile.max()
        self.profile = interp1d(t, profile/self.rate, fill_value='extrapolate', bounds_error=False)



