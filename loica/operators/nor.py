import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import least_squares
from .receiver import *

class Nor:
    color = 'orange'
    shape = 's'
    def __init__(self, input, output, alpha, a, b, K, n):
        self.alpha = alpha
        self.a = a
        self.b = b
        self.K = K
        self.n = n
        self.input = input
        self.output = output

    def __str__(self):
        return 'NOR'

    def expression_rate(self, t, dt):
        input_repressor1 = self.input[0].concentration
        input_repressor2 = self.input[1].concentration
        r1 = (input_repressor1/self.K[0])**self.n[0]
        r2 = (input_repressor2/self.K[1])**self.n[1]
        expression_rate1 = ( self.a[0] + self.b[0]*r1) / (1 + r1)
        expression_rate2 = (self.a[1] + self.b[1]*r2) / (1 + r2)
        a0 = self.alpha[0]
        a1 = expression_rate1
        a2 = expression_rate2
        a3 = self.alpha[1] * expression_rate1 * expression_rate2
        expression_rate = a1 + a2
        return expression_rate

    def forward_model(
        self,
        rep1_a=0, rep1_b=1, rep1_K=1, rep1_n=2,
        rep2_a=0, rep2_b=1, rep2_K=1, rep2_n=2,
        #alpha0=0, alpha1=0,
        a_A=1e2, b_A=0, K_A=1, n_A=2,
        a_B=1e2, b_B=0, K_B=1, n_B=2,
        Dt=0.05,
        sim_steps=10,
        A=0, B=0,
        odval=[1]*100,
        gamma=0,
        rep1_0=0, rep2_0=0, fp_0=0,
        nt=100
    ):
        fp_list,A_list,B_list,t_list = [],[],[],[]
        rep1 = rep1_0
        rep2 = rep2_0
        fp = fp_0
        for t in range(nt):
            fp_list.append(fp)
            A_list.append(A)
            B_list.append(B)
            t_list.append(t * Dt)
            od = odval[t]
            for tt in range(sim_steps):
                time = (t + tt/sim_steps) * Dt
                # Repressor 1
                a = (A/K_A)**n_A
                nextrep1 = rep1 + (od * (a_A + b_A * a) /(1 + a) - gamma*rep1) * Dt/sim_steps

                # Repressor 2
                b = (B/K_B)**n_B
                nextrep2 = rep2 + (od * (a_B + b_B * b) /(1 + b) - gamma*rep2) * Dt/sim_steps

                # Reporter output
                r1 = (rep1/od/rep1_K)**rep1_n
                r2 = (rep2/od/rep2_K)**rep2_n
                expression_rate1 = (rep1_a + rep1_b*r1) / (1 + r1)
                expression_rate2 = (rep2_a + rep2_b*r2) / (1 + r2)
                #a0 = alpha0
                a1 = expression_rate1
                a2 = expression_rate2
                #a3 = alpha1 * expression_rate1 * expression_rate2
                expression_rate = a1 + a2
                nextfp = fp + ( od * expression_rate) * Dt/sim_steps

                # Update system
                rep1,rep2,fp = nextrep1,nextrep2,nextfp

        afp = np.array(fp_list).transpose()
        AA = np.array(A_list).transpose()
        BB = np.array(B_list).transpose()
        tt = np.array(t_list).transpose()
        t = np.arange(nt) * Dt
        return afp,AA,BB,tt

    def residuals(self, df, oddf, a_A, b_A, K_A, n_A, a_B, b_B, K_B, n_B, gamma): 
        def func(x): 
            rep1_a, rep1_b, rep1_K, rep1_n = x[0:4]
            rep2_a, rep2_b, rep2_K, rep2_n = x[4:8]
            #alpha0, alpha1 = x[8:10]
            residual_list = []
            df_sorted = df.sort_values(['Sample', 'Time'])
            oddf_sorted = oddf.sort_values(['Sample', 'Time'])
            for samp_id,samp_data in df_sorted.groupby('Sample'):
                odval = oddf_sorted[oddf_sorted.Sample==samp_id].Measurement.values
                data = samp_data.Measurement.values
                p0_1 = 0
                p0_2 = data[0]
                A = samp_data.Concentration1.values[0]
                B = samp_data.Concentration2.values[0]
                t = samp_data.Time.values
                dt = np.mean(np.diff(t))
                nt = len(t)
                p,AA,BB,tt = self.forward_model(
                            rep1_a=rep1_a, rep1_b=rep1_b, rep1_K=rep1_K, rep1_n=rep1_n,
                            rep2_a=rep2_a, rep2_b=rep2_b, rep2_K=rep2_K, rep2_n=rep2_n,
                            #alpha0=alpha0, alpha1=alpha1,
                            a_A=a_A, b_A=b_A, K_A=K_A, n_A=n_A,
                            a_B=a_B, b_B=b_B, K_B=K_B, n_B=n_B,
                            Dt=dt,
                            A=A, B=B,
                            odval=odval,
                            gamma=gamma,
                            rep1_0=0, rep2_0=0, fp_0=data[0],
                            nt=nt,
                        )
                model = p.ravel()
                residual = (data[1:] - model[1:]) 
                residual_list.append(residual) 
            return np.array(residual_list).ravel()
        return func

    def characterize(self, 
            flapjack, 
            receiver1,
            receiver2, 
            nor_inverter, 
            media, 
            strain, 
            signal, 
            biomass_signal,
            gamma
            ):
        # Get biomass time series
        biomass_df = flapjack.analysis(type='Background Correct', 
                            vector=nor_inverter,
                            media=media,
                            strain=strain,
                            signal=biomass_signal,
                            biomass_signal=biomass_signal
                         )
       
        # Characterize receiver1 Hill function
        rec1 = Receiver(None, None, 0, 0, 0, 0)
        rec1.characterize(
            flapjack,
            vector=receiver1,
            media=media,
            strain=strain,
            signal=signal,
            biomass_signal=biomass_signal
        )
        self.a_A = rec1.a
        self.b_A = rec1.b
        self.K_A = rec1.K
        self.n_A = rec1.n

        # Characterize receiver2 Hill function
        rec2 = Receiver(None, None, 0, 0, 0, 0)
        rec2.characterize(
            flapjack,
            vector=receiver2,
            media=media,
            strain=strain,
            signal=signal,
            biomass_signal=biomass_signal
        )
        self.a_B = rec2.a
        self.b_B = rec2.b
        self.K_B = rec2.K
        self.n_B = rec2.n

        # Characterize inverter
        nor_inverter_df = flapjack.analysis(type='Background Correct', 
                            vector=nor_inverter,
                            media=media,
                            strain=strain,
                            signal=signal,
                            biomass_signal=biomass_signal
                         )

        # Bounds for fitting
        lower_bounds = [0]*8
        upper_bounds = [1e8, 1e8, 1e8, 8] * 2
        bounds = [lower_bounds, upper_bounds]

        '''
            rep1_a, rep1_b, rep1_K, rep1_n = x[0:4]
            rep2_a, rep2_b, rep2_K, rep2_n = x[4:8]
            alpha0, alpha1, alpha2, alpha3 = x[8:12]
        '''
        # Solve for parameters and profile
        residuals = self.residuals(
                                nor_inverter_df, biomass_df,
                                self.a_A, self.b_A, self.K_A, self.n_A,
                                self.a_B, self.b_B, self.K_B, self.n_B,
                                gamma=gamma
                            )
        res = least_squares(residuals, [1,1,1,0,1,1,1,0], bounds=bounds)
        self.res = res
        self.rep1_a, self.rep1_b, self.rep1_K, self.rep1_n = res.x[0:4]
        self.rep2_a, self.rep2_b, self.rep2_K, self.rep2_n = res.x[4:8]
        #self.alpha0, self.alpha1 = res.x[8:10]