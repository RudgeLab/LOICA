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
            od = odval[t]
            tt = t*Dt
            prof = profile[t]
            for tt in range(sim_steps):
                nextp1 = p1 + (odval[t]*profile[t] - gamma*p1) * Dt / sim_steps
                p1 = nextp1

            p1_list.append(p1)
            t_list.append([t * Dt])

        ap1 = np.array(p1_list).transpose()
        tt = np.array(t_list).transpose()
        t = np.arange(nt) * Dt
        return ap1,tt

    def residuals(self, data, p0, odval, epsilon, dt, nt, fmax): 
        def func(x): 
            gamma = x[0]
            ff = x[1::2] + x[2::2]*1j
            freqs = fftfreq(nt)
            tff = np.zeros((nt,), dtype=np.complex)
            tff[np.abs(freqs)<fmax] = ff
            profile = ifft(tff).real
            p,tt = self.forward_model(
                        Dt=dt,
                        odval=odval,
                        profile=profile,
                        nt=nt,
                        p0=p0,
                        gamma=gamma
                    )
            model = p
            residual = data - model
            tikhonov = profile
            total_variation = np.sqrt(np.abs(np.diff(profile)))
            result = np.concatenate((residual, epsilon * tikhonov))
            return result
        return func

    def characterize(self, flapjack, vector, media, strain, signal, biomass_signal, fmax, epsilon):
        expression = flapjack.analysis(media=media, 
                            strain=strain,
                            vector=vector,
                            signal=signal,
                            type='Background Correct',
                            biomass_signal=biomass_signal
                            ).sort_values(['Sample', 'Time'])
        t = expression.Time.unique()
        dt = np.diff(t).mean()
        expression = expression.groupby('Time').mean().Measurement.values

        biomass = flapjack.analysis(media=media, 
                            strain=strain,
                            vector=vector,
                            signal=biomass_signal,
                            type='Background Correct',
                            biomass_signal=biomass_signal
                            ).sort_values(['Sample', 'Time'])
        biomass = biomass.groupby('Time').mean().Measurement.values

        nt = len(biomass)
        freqs = fftfreq(nt)
        ncomps = len(freqs[np.abs(freqs)<fmax]) * 2

        # Bounds for fitting
        lower_bounds = [0] + [-1e8]*ncomps
        upper_bounds = [5] + [1e8]*ncomps
        bounds = [lower_bounds, upper_bounds]
        '''
            gamma = x[0]
            profile = x[1:]
        '''
        data = expression.ravel()
        res = least_squares(
                self.residuals(
                    data, data[0], biomass, epsilon=epsilon, dt=dt, nt=nt, fmax=fmax
                    ), 
                [0] + [1]*ncomps, 
                bounds=bounds
                )
        self.gamma = res.x[0]
        fprofile = res.x[1:]
        fcomps = fprofile[::2] + fprofile[1::2]*1j
        tff = np.zeros((nt,), dtype=np.complex)
        tff[np.abs(freqs)<fmax] = fcomps
        self.profile = interp1d(t, ifft(tff).real)



