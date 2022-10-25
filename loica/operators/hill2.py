from .operator import *
import numpy as np
from scipy.optimize import least_squares
from .receiver import *

class Hill2(Operator):
    """
    A class that represents a DNA fragment that encode a genetic operator.
    The Hill2 Operator is an abstraction of a set of two repressible or inducible promoters that
    maps an  2 inputs into an output using a Hill function.

    ...
    
    Attributes
    ----------
    input : List [Regulator | Supplement]
        The inputs of the operator that regulates the expression of the output
    output : Regulator | Reporter | List
        The output of the operator that is regulated by the input
    alpha : List
        [Basal expression rate, Regulated expression rate in MEFL/second]
    K : int | float
        Half expression input concentration in Molar 
    n : int | float
        Hill coefficient, cooperative degree (unitless)
    uri : str, optional
        SynBioHub URI
    sbol_comp : SBOL Component, optional
        SBOL Component
    name : str, optional
        Name of the operator displayed on the network representation
    color: str, optional
        Color displayed on the network representation

    Methods
    -------
    characterize(flapjack, receiver, inverter, media, strain, signal, biomass_signal, gamma)
        Parameterize the Operator model that maps Input concentration into Output expression rate
    """
    def __init__(self, input, output, alpha, K, n, name=None, uri=None, sbol_comp=None, color='orange'):
        super().__init__(output, name, uri, sbol_comp, color)
        self.alpha = alpha
        self.K = K
        self.n = n
        self.input = input

    def __str__(self):
        if self.name == None:
            return 'HILL2'
        else: return self.name

    def expression_rate(self, t, dt):
        input_repressor1 = self.input[0].concentration
        input_repressor2 = self.input[1].concentration
        r1 = (input_repressor1/self.K[0])**self.n[0]
        r2 = (input_repressor2/self.K[1])**self.n[1]
        #r12 = (input_repressor1 * input_repressor2 / self.K[0] / self.K[1])**(self.n[0] + self.n[1])
        r12 = r1 * r2
        num = self.alpha[0] + self.alpha[1] * r1 + self.alpha[2] * r2 + self.alpha[3] * r12
        denom = 1 + r1 + r2 + r12
        return num / denom

    def forward_model(
        self,
        rep1_K=1, rep1_n=2,
        rep2_K=1, rep2_n=2,
        alpha0=1, alpha1=0, alpha2=0, alpha3=0,
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
                #r12 = (rep1*rep2/od/od/rep1_K/rep2_K)**(rep1_n+rep2_n)
                r12 = r1 * r2
                num = alpha0 + alpha1 * r1 + alpha2 * r2 + alpha3 * r12
                denom = 1 + r1 + r2 + r12
                expression_rate = num / denom
                nextfp = fp + ( od * expression_rate) * Dt/sim_steps

                # Update system
                rep1,rep2,fp = nextrep1,nextrep2,nextfp

        afp = np.array(fp_list).transpose()
        AA = np.array(A_list).transpose()
        BB = np.array(B_list).transpose()
        tt = np.array(t_list).transpose()
        t = np.arange(nt) * Dt
        return afp,AA,BB,tt

    def residuals(self, df, oddf, a_A, b_A, K_A, n_A, a_B, b_B, K_B, n_B, chem1, chem2, gamma): 
        def func(x): 
            rep1_K, rep1_n = x[0:2]
            rep2_K, rep2_n = x[2:4]
            alpha0, alpha1, alpha2, alpha3 = x[4:8]
            residual_list = []
            df_sorted = df.sort_values(['Sample', 'Time'])
            oddf_sorted = oddf.sort_values(['Sample', 'Time'])
            for samp_id,samp_data in df_sorted.groupby('Sample'):
                odval = oddf_sorted[oddf_sorted.Sample==samp_id].Measurement.values
                data = samp_data.Measurement.values
                p0_1 = 0
                p0_2 = data[0]
                # Get inducer concentrations if any
                chemical1 = samp_data.Chemical_id1.values[0]
                chemical2 = samp_data.Chemical_id2.values[0]
                if chemical1 == chem1:
                    A = samp_data.Concentration1.values[0]
                elif chemical2 == chem1:
                    A = samp_data.Concentration2.values[0]
                else:
                    A = 0
                if chemical1 == chem2:
                    B = samp_data.Concentration1.values[0]
                elif chemical2 == chem2:
                    B = samp_data.Concentration2.values[0]
                else:
                    B = 0
                t = samp_data.Time.values
                dt = np.mean(np.diff(t))
                nt = len(t)
                p,AA,BB,tt = self.forward_model(
                            rep1_K=rep1_K, rep1_n=rep1_n,
                            rep2_K=rep2_K, rep2_n=rep2_n,
                            alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, alpha3=alpha3,
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
            chemical1,
            chemical2,
            nor_inverter, 
            media, 
            strain, 
            signal, 
            biomass_signal,
            gamma,
            lower_bounds=[0]*8,
            upper_bounds=[1e8, 8, 1e8, 8, 1e8, 1e8, 1e8, 1e8],
            init_x = [1,2,1,2,1,0,0,0]
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
        self.a_A = rec1.alpha[0]
        self.b_A = rec1.alpha[1]
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
        self.a_B = rec2.alpha[0]
        self.b_B = rec2.alpha[1]
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
        #lower_bounds = [0]*8
        #upper_bounds = [1e8, 8, 1e8, 8, 1e8, 1e8, 1e8, 1e8]
        bounds = [lower_bounds, upper_bounds]

        '''
            rep1_K, rep1_n = x[0:2]
            rep2_K, rep2_n = x[2:4]
            alpha0, alpha1, alpha2, alpha3 = x[4:8]
        '''
        # Solve for parameters and profile
        residuals = self.residuals(
                                nor_inverter_df, biomass_df,
                                self.a_A, self.b_A, self.K_A, self.n_A,
                                self.a_B, self.b_B, self.K_B, self.n_B,
                                chemical1, chemical2,
                                gamma=gamma
                            )
        res = least_squares(residuals, init_x, bounds=bounds)
        self.res = res
        self.rep1_K, self.rep1_n = res.x[0:2]
        self.rep2_K, self.rep2_n = res.x[2:4]
        self.alpha0, self.alpha1, self.alpha2, self.alpha3 = res.x[4:8]
