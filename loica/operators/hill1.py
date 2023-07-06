from .operator import *
import numpy as np
from scipy.optimize import least_squares
from scipy.integrate import solve_ivp, odeint
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from flapjack import *
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
        n_i,
        K_i,
        od,
        profile,
        rec_profile,
        gamma,
        p0_1, p0_2,
        tt,
        sim_steps=10
    ):
        # od, rec_profile = callable as function of time
        def dydt(y, t):
            #a = (A/K_A)**n_A
            #p = (p1/od(t)/K_i)**n_i
            #dp1dt = od(t) * (a_A + b_A * a) /(1 + a) - gamma*p1
            rec = max(rec_profile(t), 1)
            odval = max(od(t), 0.001)
            p1 = rec / odval
            p = (p1 / K_i)**n_i
            dp2dt = odval * profile(t) * (a_j + b_j*p) / ( 1 + p )
            return dp2dt
        #res = solve_ivp(fun=dydt, t_span=[0, np.max(tt)], y0=np.array([p0_1, p0_2], dtype=float), t_eval = tt)
        #y = res.sol(tt)
        y = odeint(dydt, p0_2, tt, rtol=1e-3)
        p2 = y[:,0]

        '''
        p2_list,A_list,t_list = [],[],[]
        p2 = p0_2
        p2_list.append(p2)
        t_list.append(tt[0])
        nt = len(tt)
        for it in range(1, nt):
            t = tt[it]
            Dt = t - tt[it-1]
            if od(t)>0:
                p2_list.append(p2)
                t_list.append(t)
                for st in range(sim_steps):
                    time = (t + st/sim_steps) * Dt
                    rec = max(rec_profile(time), 1)            
                    p1 = rec / od(time)
                    p = (p1/K_i)**n_i
                    nextp2 = p2 + ( od(time) * (a_j + b_j*p) / ( 1 + p )) * Dt/sim_steps
                    p2 = nextp2

        p2 = np.array(p2_list).transpose()
        AA = np.array(A_list).transpose()
        tt = np.array(t_list).transpose()
        t = np.arange(nt) * Dt
        '''
        return p2,tt

    def residuals(self, df, oddf, rec_df, profile, gamma, plotting=False): 
        def func(x): 
            self.data = []
            self.model = []
            b_j = x[3] #np.exp(x[3])
            a_j = x[2] #b_j / np.exp(x[2]) 
            n_i = x[0] #np.exp(x[0])
            K_i = x[1] #np.exp(x[1])  
            print(a_j, b_j, K_i, n_i)
            
            residual_array = np.array((0,))
            
            df_sorted = df.sort_values(['Sample', 'Time'])
            oddf_sorted = oddf.sort_values(['Sample', 'Time'])
            df_sorted = df_sorted[oddf_sorted.Measurement>0]
            oddf_sorted = oddf_sorted[oddf_sorted.Measurement>0]

            for samp_id,samp_data in df_sorted.groupby('Sample'):
                odval = oddf_sorted[oddf_sorted.Sample==samp_id].Measurement.values.astype(float)
                tt = oddf_sorted[oddf_sorted.Sample==samp_id].Time.values.astype(float)
                od = interp1d(tt, odval, bounds_error=False, fill_value='extrapolate')
                p0_1 = 0
                p0_2 = samp_data.Measurement.values[0]
                A = samp_data.Concentration1.values[0]
                #if np.isnan(A):
                #    A = 0
                #    rec_profile = rec_df[rec_df.Concentration1.isna()].sort_values('Time').groupby('Time').mean().Measurement.values
                #    rec_profile_t = rec_df[rec_df.Concentration1.isna()].sort_values('Time').groupby('Time').mean().index
                #else:
                if ~np.isnan(A):
                    rec_profile = rec_df[rec_df.Concentration1==A].sort_values('Time').groupby('Time').mean().Measurement.values
                    rec_profile_t = rec_df[rec_df.Concentration1==A].sort_values('Time').groupby('Time').mean().index
                    #print(A)
                    #print(rec_profile_t)
                    #print(rec_profile)
                    rec_profile = interp1d(rec_profile_t[rec_profile>0], rec_profile[rec_profile>0], bounds_error=False, fill_value='extrapolate')
                    p,tt = self.forward_model(
                                a_j=a_j,
                                b_j=b_j,
                                n_i=n_i,
                                K_i=K_i,
                                od=od,
                                profile=profile,
                                rec_profile=rec_profile,
                                gamma=gamma,
                                p0_1=p0_1, p0_2=p0_2,
                                tt=tt
                            )
                    data = samp_data.Measurement.values
                    self.data.append(data)
                    model = p.ravel()
                    self.model.append(model)

                    if plotting:
                        plt.plot(tt, model, '--')
                        plt.plot(tt, data)
                    residual_array = np.append(residual_array, data[1:] - model[1:])
            print(np.sum(residual_array * residual_array))
            return residual_array.ravel()
        return func

    def characterize(self, 
            flapjack, 
            receiver, 
            inverter, 
            study,
            media, 
            strain, 
            signal, 
            biomass_signal,
            gamma
            ):
        # Get biomass time series
        biomass_df = flapjack.analysis(type='Background Correct', 
                            study=study,
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
        
        self.a_A = receiver.alpha[0]
        self.b_A = receiver.alpha[1]
        self.K_A = receiver.K
        self.n_A = receiver.n
        '''
        
        rec_df = flapjack.analysis(vector=receiver,
                            study=study,
                            signal=signal,
                            media=media,
                            strain=strain,
                            type='Background Correct',
                            biomass_signal=biomass_signal)

        inv_df = flapjack.analysis(vector=inverter,
                            study=study,
                            signal=signal,
                            media=media,
                            strain=strain,
                            type='Background Correct',
                            biomass_signal=biomass_signal)
        inv_er_df = expression_rate_inverse(inv_df, biomass_df)
        inv_er_df = inv_er_df[inv_er_df.Concentration1.isna()]
        inv_er_val = inv_er_df.groupby('Time').mean().Rate.values
        inv_er_t = inv_er_df.groupby('Time').mean().index
        
        profile = interp1d(inv_er_t, inv_er_val/inv_er_val.max(), bounds_error=False, fill_value='extrapolate')


        # Characterize inverter
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
        initx = [n, K, a, b] 
        lower_bounds = [1, 0, 0, 0]
        upper_bounds = [6, 1e8, 1e8, 1e8]
        bounds = [lower_bounds, upper_bounds]

        func = self.residuals(inv_df,
                        biomass_df,
                        rec_df,
                        profile,
                        gamma=gamma,
                        plotting=False)
        # Solve for parameters
        res = least_squares(fun=func, 
                            x0=initx,
                            bounds=bounds,
                            diff_step=[0.1,0.1,0.1,0.1]
                            #ftol=1e-3
        )

        func = self.residuals(inv_df,
                        biomass_df,
                        rec_df,
                        profile,
                        gamma=gamma,
                        plotting=True)
        func(res.x)

        plt.show()

        print(res)
        self.res = res
        self.n = np.exp(res.x[0])
        self.K = np.exp(res.x[1])
        self.alpha[1] = np.exp(res.x[3])
        self.alpha[0] = self.alpha[1] / np.exp(res.x[2])
        self.profile = profile
      
