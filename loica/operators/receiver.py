import numpy as np
from numpy.fft import fft, ifft, fftfreq
from scipy.optimize import least_squares
from scipy.interpolate import interp1d

from .source import *
from flapjack import *

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
            for tt in range(sim_steps):
                a = (A/K_A)**n_A
                prof = profile( (t + tt / sim_steps) * Dt)
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

    def residuals(self, data, p0, profile, A, odval, epsilon, dt, nt, fmax): 
        def func(x): 
            K_A = x[0]
            n_A = x[1]        
            
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
            return residual
        return func

    def characterize(self, flapjack, vector, media, strain, signal, biomass_signal, analyte, fmax, epsilon):
        source = Source(None, 0)
        source.characterize(
                    flapjack,
                    vector=vector,
                    media=media,
                    strain=strain,
                    signal=signal,
                    biomass_signal=biomass_signal,
                    fmax=fmax,
                    epsilon=epsilon
                )
        self.rate = source.rate
        self.profile = source.profile

        curve = flapjack.analysis(media=media, 
                            strain=strain,
                            vector=vector,
                            signal=signal,
                            function='Mean Expression',
                            type='Induction Curve',
                            analyte=analyte,
                            biomass_signal=biomass_signal
                            )
        params,std = fit_curve(hill, curve, x='Concentration', y='Expression')
        self.params = params
        self.std = std
        self.curve = curve
        
