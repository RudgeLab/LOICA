import numpy as np
from numpy.fft import fft, ifft, fftfreq
from scipy.optimize import least_squares
from scipy.interpolate import interp1d

class Receiver:
    def __init__(self, inducer, output, a, b, K, n, profile=None):
        if profile:
            self.profile = profile
        else:
            def profile(t):
                return 1
            self.profile = profile
        self.a = a
        self.b = b
        self.K = K
        self.n = n
        self.receptor = 0
        self.inducer = inducer
        self.output = output
        
    def expression_rate(self, t, dt):
        inducer = self.inducer.concentration
        i = (inducer/self.K)**self.n
        expression_rate = self.profile(t) * ( self.a + self.b*i ) / (1 + i)
        return expression_rate

    def forward_model(
            self,
            K_A=1,
            n_A=2,
            Dt=0.05,
            sim_steps=10,
            A=[0],
            odval=[1]*100,
            profile=[1]*100,
            gamma=0,
            p0=0,
            nt=100
        ):
        p1_list,A_list,t_list = [],[],[]
        p1 = np.zeros_like(A) + p0
        for t in range(nt):
            od = odval[t]
            tt = t*Dt
            prof = profile[t]
            for tt in range(sim_steps):
                a = (A/K_A)**n_A
                nextp1 = p1 + (od * prof * a/(1 + a) - gamma*p1) * Dt/sim_steps
                p1 = nextp1

            p1_list.append(p1)
            A_list.append(A)
            t_list.append([t * Dt]*len(A))

        ap1 = np.array(p1_list).transpose()
        AA = np.array(A_list).transpose()
        tt = np.array(t_list).transpose()
        t = np.arange(nt) * Dt
        return ap1,AA,tt

    def residuals(self, data, p0, A, odval, epsilon, dt, nt, fmax): 
        def func(x): 
            K_A = x[0]
            n_A = x[1]        
            ff = x[2::2] + x[3::2]*1j
            freqs = fftfreq(nt)
            tff = np.zeros((nt,), dtype=np.complex)
            tff[np.abs(freqs)<fmax] = ff
            profile = ifft(tff).real
            
            p,AA,tt = self.forward_model(
                        K_A=K_A,
                        n_A=n_A,
                        Dt=dt,
                        A=A,
                        odval=odval,
                        profile=profile,
                        nt=nt,
                        p0=p0
                    )
            model = p.ravel()
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
                            )
        # Time points and interval
        t = np.linspace(0, 24, 100) #expression.Time.unique()
        print(t)
        dt = np.diff(t).mean()
        # Inducer concentrations
        A = expression.groupby('Concentration1').mean().index.values
        # Group and average data
        expression = expression.sort_values(['Sample', 'Concentration1', 'Time'])
        expression = expression.groupby(['Concentration1', 'Time']).mean().Measurement.values

        biomass = flapjack.analysis(media=media, 
                            strain=strain,
                            vector=vector,
                            signal=biomass_signal,
                            type='Background Correct',
                            biomass_signal=biomass_signal
                            )
        biomass = biomass.sort_values(['Sample', 'Concentration1', 'Time'])        
        biomass = biomass.groupby('Time').mean().Measurement.values

        nt = len(t)
        freqs = fftfreq(nt)
        ncomps = len(freqs[np.abs(freqs)<fmax]) * 2

        # Bounds for fitting
        lower_bounds = [0]*2 + [-1e8]*ncomps
        upper_bounds = [1e2, 4] + [1e8]*ncomps
        bounds = [lower_bounds, upper_bounds]
        '''
            K_A = x[0]
            n_A = x[1]
            profile = x[2:]
        '''

        data = expression.ravel()
        res = least_squares(
                self.residuals(
                    data, data[0], A, biomass, epsilon=epsilon, dt=dt, nt=nt, fmax=fmax
                    ), 
                [0, 0] + [1]*ncomps, 
                bounds=bounds
                )
        self.K = res.x[0]
        self.n = res.x[1]
        fprofile = res.x[2:]
        fcomps = fprofile[::2] + fprofile[1::2]*1j
        tff = np.zeros((nt,), dtype=np.complex)
        tff[np.abs(freqs)<fmax] = fcomps
        self.profile = interp1d(t, ifft(tff).real)

