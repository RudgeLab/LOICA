import numpy as np
from scipy.interpolate import interp1d
from numpy.fft import fft, ifft, fftfreq
from scipy.optimize import least_squares
from .receiver import *

class Not:
    def __init__(self, input, output, a, b, K, n, profile=None):
        self.a = a
        self.b = b
        self.K = K
        self.n = n
        self.input = input
        self.output = output
        if not profile:
            def profile(t):
                return 1
        self.profile = profile
        
    def expression_rate(self, t, dt):
        input_repressor = self.input.concentration
        r = (input_repressor/self.K)**self.n
        expression_rate = self.profile(t) * ( self.a + self.b*r ) / (1 + r)
        return expression_rate

    def forward_model_growth(
        self,
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

    def forward_model(
        self,
        b_j,
        n_i=2,
        K_i=1,
        alpha_A=1e2,
        K_A=1,
        n_A=2,
        Dt=0.05,
        sim_steps=10,
        A=[0],
        muval=[1]*100, odval=[1]*100,
        profile_j=[1]*100,
        profile_A=[1]*100,
        gamma=0,
        p0_1=0, p0_2=0,
        nt=100
    ):
        p2_list,A_list,t_list = [],[],[]
        p1 = np.zeros_like(A) + p0_1
        p2 = np.zeros_like(A) + p0_2
        for t in range(nt):
            p2_list.append(p2)
            A_list.append(A)
            t_list.append([t * Dt]*len(A))
            od = odval[t]
            mu = muval[t]
            prof_A = profile_A(t * Dt)
            prof_j = profile_j[t]
            for tt in range(sim_steps):
                time = (t + tt/sim_steps) * Dt
                a = (A/K_A)**n_A
                nextp1 = p1 + (prof_A * a/(1 + a) - gamma*p1 - mu*p1) * Dt/sim_steps
                p = (p1/K_i)**n_i
                nextp2 = p2 + ( od * prof_j * (1 + b_j*p) / ( 1 + p )) * Dt/sim_steps
                p1,p2 = nextp1,nextp2


        ap2 = np.array(p2_list).transpose()
        AA = np.array(A_list).transpose()
        tt = np.array(t_list).transpose()
        t = np.arange(nt) * Dt
        return ap2,AA,tt

    def residuals_growth(self, data, od0, epsilon, dt, nt, fmax): 
        def func(x): 
            ff = x[::2] + x[1::2]*1j
            freqs = fftfreq(nt)
            tff = np.zeros((nt,), dtype=np.complex)
            tff[np.abs(freqs)<fmax] = ff
            muval = ifft(tff).real

            od,tt = self.forward_model_growth(
                        Dt=dt,
                        muval=muval,
                        od0=od0,
                        nt=nt
                    )
            model = np.tile(od, int(len(data)/nt))
            residual = data - model
            tikhonov = muval
            total_variation = np.sqrt(np.abs(np.diff(muval)))
            result = np.concatenate((residual, epsilon * tikhonov))
            return residual
        return func

    def residuals(self, data, p0_1, p0_2, profile_A, K_A, n_A, A, muval, odval, epsilon, Dt, nt, fmax): 
        def func(x): 
            n_i = x[0]
            K_i = x[1]
            gamma = x[2]
            b_j = x[3]
            ff = x[4::2] + x[5::2]*1j
            freqs = fftfreq(nt)
            tff = np.zeros((nt,), dtype=np.complex)
            tff[np.abs(freqs)<fmax] = ff
            profile_j = ifft(tff).real
            
            p,AA,tt = self.forward_model(
                        b_j,
                        n_i=n_i,
                        K_i=K_i,
                        K_A=K_A,
                        n_A=n_A,
                        Dt=Dt,
                        A=A,
                        muval=muval, odval=odval,
                        profile_j=profile_j,
                        profile_A=profile_A,
                        gamma=gamma,
                        nt=nt,
                        p0_1=p0_1, p0_2=p0_2
                    )
            model = p.ravel()
            residual = data - model
            tikhonov = profile_j
            result = np.concatenate((residual, epsilon * tikhonov))
            return residual
        return func

    def characterize_growth(self,
            flapjack, 
            inverter, 
            media, 
            strain, 
            biomass_signal, 
            fmax, 
            epsilon
            ):
        # Characterize growth rate profile
        biomass = flapjack.analysis(type='Background Correct', 
                            vector=inverter,
                            media=media,
                            strain=strain,
                            signal=biomass_signal,
                            biomass_signal=biomass_signal
                         )
        biomass = biomass.sort_values('Time')
        t = biomass.Time.unique()
        dt = np.mean(np.diff(t))
        nt = len(t)
        freqs = fftfreq(nt)
        ncomps = len(freqs[np.abs(freqs)<fmax]) * 2

        biomass = biomass.groupby('Time').mean().Measurement.values
        lower_bounds = [-50]*ncomps
        upper_bounds = [50]*ncomps
        bounds = [lower_bounds, upper_bounds]

        data = biomass
        res = least_squares(
                self.residuals_growth(data, 0.01, epsilon=epsilon, dt=dt, nt=nt, fmax=fmax), 
                [1]*ncomps, 
                bounds=bounds
                )
        fprofile = res.x
        fcomps = fprofile[::2] + fprofile[1::2]*1j
        freqs = fftfreq(nt)
        tff = np.zeros((nt,), dtype=np.complex)
        tff[np.abs(freqs)<fmax] = fcomps

        self.mu_profile = ifft(tff).real
        self.biomass = biomass


    def characterize(self, 
            flapjack, 
            receiver, 
            inverter, 
            media, 
            strain, 
            signal, 
            biomass_signal, 
            analyte, 
            fmax, 
            epsilon
            ):
        # Get growth rate profile
        self.characterize_growth(flapjack, 
            inverter, 
            media, 
            strain, 
            biomass_signal, 
            fmax, 
            epsilon)

        # Characterize receiver profile and Hill function
        rec = Receiver(None, None, 0, 0, 0, 0)
        rec.characterize(
            flapjack,
            vector=receiver,
            media=media,
            strain=strain,
            signal=signal,
            biomass_signal=biomass_signal,
            analyte=analyte,
            fmax=0.1,
            epsilon=0
        )
        alpha_A = rec.a
        profile_A = rec.profile
        K_A = rec.K
        n_A = rec.n

        # Characterize inverter
        inverter = flapjack.analysis(type='Background Correct', 
                            vector=inverter,
                            media=media,
                            strain=strain,
                            signal=signal,
                            biomass_signal=biomass_signal
                         )
        inverter = inverter.sort_values(['Sample', 'Concentration1', 'Time'])
        t = inverter.Time.unique()
        dt = np.mean(np.diff(t))
        nt = len(t)
        # Parameters for Fourier basis
        freqs = fftfreq(nt)
        ncomps = len(freqs[np.abs(freqs)<fmax]) * 2

        A = inverter.groupby('Concentration1').mean().index.values
        inverter = inverter.groupby(['Concentration1', 'Time']).mean().Measurement.values

        # Bounds for fitting
        lower_bounds = [0]*4 + [-1e10]*ncomps
        upper_bounds = [4, 1e8, 3, 1e-3]  + [1e10]*ncomps
        bounds = [lower_bounds, upper_bounds]

        '''
            n_i = x[0]
            K_i = x[1]
            gamma = x[2]
            b_j = x[3]
            profile_j = x[4:]
        '''
        # Solve for parameters and profile
        data = inverter.ravel()
        residuals = self.residuals(
                                data, 
                                0, data[0],
                                profile_A, 
                                K_A, n_A, A,
                                self.mu_profile, self.biomass, 
                                epsilon=0, Dt=dt, nt=nt,
                                fmax=fmax
                            )
        res = least_squares(residuals, [1,1,1,0] + [1]*ncomps, bounds=bounds)
        self.n = res.x[0]
        self.K = res.x[1]
        self.gamma = res.x[2]
        self.b = res.x[3]
        fprofile = res.x[4:]
        fcomps = fprofile[::2] + fprofile[1::2]*1j
        freqs = fftfreq(nt)
        tff = np.zeros((nt,), dtype=np.complex)
        tff[np.abs(freqs)<fmax] = fcomps
        profile = ifft(tff).real
        self.a = profile.max()
        self.profile = interp1d(t, profile/self.a, fill_value='extrapolate', bounds_error=False)
