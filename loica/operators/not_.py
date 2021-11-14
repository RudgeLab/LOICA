import numpy as np
from scipy.interpolate import interp1d
from numpy.fft import fft, ifft, fftfreq
from scipy.optimize import least_squares
from .receiver import *

class Not:
    color = 'skyblue'
    shape = 's'
    def __init__(self, input, output, a, b, K, n, uri=None, sbol_comp=None):
        self.a = a
        self.b = b
        self.K = K
        self.n = n
        self.input = input
        self.output = output
        self.uri = uri
        self.sbol_comp = sbol_comp

    def __str__(self):
        return 'NOT'
        
    def expression_rate(self, t, dt):
        input_repressor = self.input.concentration
        r = (input_repressor/self.K)**self.n
        expression_rate = ( self.a + self.b*r ) / (1 + r)
        return expression_rate

    def forward_model(
        self,
        a_j,
        b_j,
        n_i=2,
        K_i=1,
        a_A=1e2,
        b_A=0,
        K_A=1,
        n_A=2,
        Dt=0.05,
        sim_steps=10,
        A=[0],
        odval=[1]*100,
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
            for tt in range(sim_steps):
                time = (t + tt/sim_steps) * Dt
                a = (A/K_A)**n_A
                nextp1 = p1 + (od * (a_A + b_A * a) /(1 + a) - gamma*p1) * Dt/sim_steps
                p = (p1/od/K_i)**n_i
                nextp2 = p2 + ( od * (a_j + b_j*p) / ( 1 + p )) * Dt/sim_steps
                p1,p2 = nextp1,nextp2

        ap2 = np.array(p2_list).transpose()
        AA = np.array(A_list).transpose()
        tt = np.array(t_list).transpose()
        t = np.arange(nt) * Dt
        return ap2,AA,tt

    def residuals(self, data, p0_1, p0_2, a_A, b_A, K_A, n_A, A, odval, gamma, Dt, t): 
        def func(x): 
            nt = len(t)
            n_i = x[0]
            K_i = x[1]  
            a_j = x[2]          
            b_j = x[3]         
            p,AA,tt = self.forward_model(
                        a_j=a_j,
                        b_j=b_j,
                        n_i=n_i,
                        K_i=K_i,
                        a_A=a_A,
                        b_A=b_A,
                        K_A=K_A,
                        n_A=n_A,
                        Dt=Dt,
                        A=A,
                        odval=odval,
                        gamma=gamma,
                        nt=nt,
                        p0_1=p0_1, p0_2=p0_2
                    )
            model = p.ravel()
            residual = (data[1:] - model[1:]) 
            return residual
        return func

    def characterize(self, 
            flapjack, 
            receiver, 
            inverter, 
            media, 
            strain, 
            signal, 
            biomass_signal,
            gamma
            ):
        # Get biomass time series
        biomass = flapjack.analysis(type='Background Correct', 
                            vector=inverter,
                            media=media,
                            strain=strain,
                            signal=biomass_signal,
                            biomass_signal=biomass_signal
                         )
        biomass = biomass.sort_values('Time')
        biomass = biomass.groupby('Time').mean().Measurement.values
        
        # Characterize receiver profile and Hill function
        rec = Receiver(None, None, 0, 0, 0, 0)
        rec.characterize(
            flapjack,
            vector=receiver,
            media=media,
            strain=strain,
            signal=signal,
            biomass_signal=biomass_signal
        )
        self.a_A = rec.a
        self.b_A = rec.b
        self.K_A = rec.K
        self.n_A = rec.n

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

        A = inverter.groupby('Concentration1').mean().index.values
        inverter = inverter.groupby(['Concentration1', 'Time']).mean().Measurement.values

        # Bounds for fitting
        lower_bounds = [0]*4
        upper_bounds = [8, 1e8, 1e8, 1e8]
        bounds = [lower_bounds, upper_bounds]

        '''
            n_i = x[0]
            K_i = x[1]
            a_j = x[2]
            b_j = x[3]
        '''
        # Solve for parameters and profile
        data = inverter.ravel()
        residuals = self.residuals(
                                data, 
                                0, data[0], 
                                self.a_A, self.b_A, self.K_A, self.n_A, A,
                                biomass, 
                                gamma=gamma,
                                Dt=dt, t=t
                            )
        res = least_squares(residuals, [1,1,1,0], bounds=bounds)
        self.res = res
        self.n = res.x[0]
        self.K = res.x[1]
        self.a = res.x[2]
        self.b = res.x[3]