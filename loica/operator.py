import numpy as np
from numpy.fft import fft, ifft, fftfreq
from scipy.optimize import least_squares

class Operator:
    def expression_rate(self):
        return 0


class Source(Operator):
    def __init__(self, output, rate):
        self.rate = rate
        self.output = output

    def expression_rate(self):
        return self.rate

    def simulate(
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
            p,tt = self.simulate(
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
        dt = np.diff(expression.Time.unique()).mean()
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
        self.profile = ifft(tff).real

class Receiver(Operator):
    def __init__(self, receptor, inducer, output, a, b, K, n):
        self.a = a
        self.b = b
        self.K = K
        self.n = n
        self.receptor = 0
        self.inducer = inducer
        self.receptor = receptor
        self.output = output
        
    def expression_rate(self):
        receptor = self.receptor.concentration
        inducer = self.inducer.concentration
        complex = receptor * inducer
        i = (complex/self.K)**self.n
        expression_rate = ( self.a + self.b*i ) / (1 + i)
        return expression_rate

class Not(Operator):
    def __init__(self, input, output, a, b, K, n):
        self.a = a
        self.b = b
        self.K = K
        self.n = n
        self.input = input
        self.output = output
        
    def expression_rate(self):
        input_repressor = self.input.concentration
        r = (input_repressor/self.K)**self.n
        expression_rate = ( self.a + self.b*r ) / (1 + r)
        return expression_rate


class Buffer(Operator):
    def __init__(self, input, output, a, b, K, n):
        self.a = a
        self.b = b
        self.K = K
        self.n = n
        self.input = input
        self.output = output
        
    def expression_rate(self):
        input_inducer = self.input.concentration
        i = (input_inducer/self.K)**self.n
        expression_rate = ( self.a + self.b*i ) / (1 + i)
        return expression_rate


