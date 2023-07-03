from .operator import *
import numpy as np
from scipy.optimize import least_squares
from .receiver import *

class Hill1(Operator):
    """
    A class that represents a DNA fragment that encode a genetic operator.
    The Hill1 Operator is an abstraction of a repressible or inducible promoter that
    maps an input into an output using a Hill function.

    ...
    
    Attributes
    ----------
    input : Regulator | Supplement
        The input of the operator that regulates the expression of the output
    output : Regulator | Reporter
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

    def __init__(self, input, output, alpha, K, n, name=None, uri=None, sbol_comp=None, color='skyblue'):
        super().__init__(output, name, uri, sbol_comp, color)
        self.alpha = alpha
        self.K = K
        self.n = n
        self.input = input

    def __str__(self):
        if self.name == None:
            return 'HILL1'
        else: return self.name
        
    def expression_rate(self, t, dt):
        input_repressor = self.input.concentration
        r = (input_repressor/self.K)**self.n
        expression_rate = ( self.alpha[0] + self.alpha[1]*r ) / (1 + r)
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
        A=0,
        odval=[1]*100,
        gamma=0,
        p0_1=0, p0_2=0,
        nt=100
    ):
        p2_list,A_list,t_list = [],[],[]
        p1 = p0_1
        p2 = p0_2
        for t in range(nt):
            p2_list.append(p2)
            A_list.append(A)
            t_list.append(t * Dt)
            od = odval[t]
            for tt in range(sim_steps):
                time = (t + tt/sim_steps) * Dt
                a = (A/K_A)**n_A
                nextp1 = p1 + ((a_A + b_A * a) /(1 + a) - gamma*p1) * Dt/sim_steps
                p = (p1/K_i)**n_i
                nextp2 = p2 + ( od * (a_j + b_j*p) / ( 1 + p )) * Dt/sim_steps
                p1,p2 = nextp1,nextp2

        ap2 = np.array(p2_list).transpose()
        AA = np.array(A_list).transpose()
        tt = np.array(t_list).transpose()
        t = np.arange(nt) * Dt
        return ap2,AA,tt

    def residuals(self, df, oddf, a_A, b_A, K_A, n_A, gamma): 
        def func(x): 
            b_j = np.exp(x[3])
            a_j = b_j/ np.exp(x[2]) 
            n_i = np.exp(x[0])
            K_i = np.exp(x[1])  
            print(a_j, b_j, K_i, n_i)
            residual_list = []
            df_sorted = df.sort_values(['Sample', 'Time'])
            oddf_sorted = oddf.sort_values(['Sample', 'Time'])
            for samp_id,samp_data in df_sorted.groupby('Sample'):
                odval = oddf_sorted[oddf_sorted.Sample==samp_id].Measurement.values
                data = samp_data.Measurement.values
                p0_1 = 0
                p0_2 = data[0]
                A = samp_data.Concentration1.values[0]
                if np.isnan(A):
                    A = 0
                t = samp_data.Time.values
                dt = np.mean(np.diff(t))
                nt = len(t)
                p,AA,tt = self.forward_model(
                            a_j=a_j,
                            b_j=b_j,
                            n_i=n_i,
                            K_i=K_i,
                            a_A=a_A,
                            b_A=b_A,
                            K_A=K_A,
                            n_A=n_A,
                            Dt=dt,
                            A=A,
                            odval=odval,
                            gamma=gamma,
                            nt=nt,
                            p0_1=p0_1, p0_2=p0_2
                        )
                model = p.ravel()
                residual = (data[1:] - model[1:]) 
                residual_list.append(residual) 
            residual_array = np.array(residual_list).ravel()
            return residual_array
        return func

    def characterize(self, 
            flapjack, 
            receiver: Operator, 
            inverter, 
            media, 
            strain, 
            signal, 
            biomass_signal,
            gamma
            ):
        # Get biomass time series
        biomass_df = flapjack.analysis(type='Background Correct', 
                            vector=inverter,
                            media=media,
                            strain=strain,
                            signal=biomass_signal,
                            biomass_signal=biomass_signal
                         )
        
        # Characterize receiver profile and Hill function
        '''
        rec = Receiver(None, None, 0, 0, 0, 0)
        rec.characterize(
            flapjack,
            vector=receiver,
            media=media,
            strain=strain,
            signal=signal,
            biomass_signal=biomass_signal
        )
        '''
        self.a_A = receiver.alpha[0]
        self.b_A = receiver.alpha[1]
        self.K_A = receiver.K
        self.n_A = receiver.n

        # Characterize inverter
        inverter_df = flapjack.analysis(type='Background Correct', 
                            vector=inverter,
                            media=media,
                            strain=strain,
                            signal=signal,
                            biomass_signal=biomass_signal
                         )
        '''
            b_j = np.exp(x[3])
            a_j = b_j/ np.exp(x[2]) 
            n_i = np.exp(x[0])
            K_i = np.exp(x[1])
        '''
        a = self.alpha[0]
        b = self.alpha[1]
        K = self.K
        n = self.n
        initx = [np.log(n), np.log(K), np.log(b/a), np.log(b)] 

        # Solve for parameters
        print(inverter_df)
        print(biomass_df)
        res = least_squares(self.residuals(
                                inverter_df,
                                biomass_df, 
                                self.a_A, self.b_A, self.K_A, self.n_A,
                                gamma=gamma
                            ), 
                            initx,
        )

        self.res = res
        self.n = np.exp(res.x[0])
        self.K = np.exp(res.x[1])
        self.alpha[1] = np.exp(res.x[3])
        self.alpha[0] = self.alpha[1] / np.exp(res.x[2])
      
