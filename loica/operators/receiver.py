import numpy as np
from numpy.fft import fft, ifft, fftfreq
from scipy.optimize import least_squares
from scipy.interpolate import interp1d

from .source import *
from flapjack import *

class Receiver:
    color = 'orange'
    def __init__(self, input, output, a, b, K, n):
        self.a = a
        self.b = b
        self.K = K
        self.n = n
        self.input = input
        self.output = output

    def __str__(self):
        return 'REC'
        
    def expression_rate(self, t, dt):
        inducer = self.input.concentration
        i = (inducer/self.K)**self.n
        expression_rate = ( self.a + self.b*i ) / (1 + i)
        return expression_rate

    def forward_model(
            self,
            a=0,
            b=1,
            K_A=1,
            n_A=2,
            Dt=0.05,
            sim_steps=10,
            A=[0],
            odval=[1]*100,
            gamma=0,
            p0=0,
            nt=100
        ):
        p1_list,A_list,t_list = [],[],[]
        p1 = p0
        for t in range(nt):
            p1_list.append(p1)
            A_list.append(A)
            t_list.append(t * Dt)
            od = odval[t]
            tt = t*Dt
            for tt in range(sim_steps):
                aa = (A/K_A)**n_A
                nextp1 = p1 + (od * (a + b*aa)/(1 + aa) - gamma*p1) * Dt/sim_steps
                p1 = nextp1

        ap1 = np.array(p1_list).transpose()
        AA = np.array(A_list).transpose()
        tt = np.array(t_list).transpose()
        t = np.arange(nt) * Dt
        return ap1,AA,tt

    def residuals(self, df, oddf): 
        def func(x): 
            a = x[0]
            b = x[1]        
            K_A = x[2]
            n_A = x[3]        
            residual_list = []
            df_sorted = df.sort_values(['Sample', 'Time'])
            oddf_sorted = oddf.sort_values(['Sample', 'Time'])
            for samp_id,samp_data in df_sorted.groupby('Sample'):
                odval = oddf_sorted[oddf_sorted.Sample==samp_id].Measurement.values
                data = samp_data.Measurement.values
                p0 = data[0]
                A = samp_data.Concentration1.values[0]
                t = samp_data.Time.values
                dt = np.mean(np.diff(t))
                nt = len(t)
                p,AA,tt = self.forward_model(
                            a=a,
                            b=b,
                            K_A=K_A,
                            n_A=n_A,
                            Dt=dt,
                            A=A,
                            odval=odval,
                            nt=nt,
                            p0=p0
                        )
                model = p.ravel()
                residual = data - model
                residual_list.append(residual) 
            return np.array(residual_list).ravel()
        return func


    def characterize(self, flapjack, vector, media, strain, signal, biomass_signal):
        expression_df = flapjack.analysis(media=media, 
                            strain=strain,
                            vector=vector,
                            signal=signal,
                            type='Background Correct',
                            biomass_signal=biomass_signal
                            )

        biomass_df = flapjack.analysis(media=media, 
                            strain=strain,
                            vector=vector,
                            signal=biomass_signal,
                            type='Background Correct',
                            biomass_signal=biomass_signal
                            )

        # Bounds for fitting (a,b,K,n)
        lower_bounds = [0, 0, 0, 0]
        upper_bounds = [1e8, 1e8, 1e8, 6]
        bounds = [lower_bounds, upper_bounds]
        '''
            a = x[0]
            b = x[1]
            K_A = x[2]
            n_A = x[3]
        '''

        res = least_squares(
                self.residuals(
                    expression_df, biomass_df
                    ), 
                [0, 1, 1, 1], 
                bounds=bounds
                )
        self.res = res
        self.a = res.x[0]
        self.b = res.x[1]
        self.K = res.x[2]
        self.n = res.x[3]