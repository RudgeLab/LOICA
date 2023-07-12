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

    def __init__(self, input, output, alpha, K, n, profile=None, name=None, uri=None, sbol_comp=None, color='skyblue'):
        super().__init__(output, name, uri, sbol_comp, color)
        if not profile:
           def profile(t):
            return 1
        self.profile = profile
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
        expression_rate = self.profile(t) * ( self.alpha[0] + self.alpha[1]*r ) / (1 + r)
        return expression_rate

    def make_profile(self, heights, t, normalize=False):
        n_gaussians = len(heights)
        means = np.linspace(t.min(), t.max(), n_gaussians)
        var = (t.max()-t.min())/n_gaussians
        profile = np.zeros_like(t)
        for mean,height in zip(means, heights):
            gaussian = height * np.exp(-(t-mean)*(t-mean) / var / 2) / np.sqrt(2 * np.pi * var)
            profile = profile + gaussian
        if normalize:
            profile = profile / profile.max()
        profile = interp1d(t, profile, bounds_error=False, fill_value='extrapolate')
        return profile


    def forward_model(
        self,
        a, b, K, n,
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
            p2 = y
            rec = max(rec_profile(t), 0)
            odval = max(od(t), 0.001)
            p1 = rec / odval
            p = (p1 / K)**n
            dp2dt = odval * profile(t) * (a + b * p) / ( 1 + p )
            #sigma = profile(t) / (a_0 * k_sigma  - k_sigma * profile(t))
            #dp2dt = odval * (a_0 * profile(t) + a_1 * k_1 * profile(t) * p1**n) / (1 + profile(t) + k_1 * profile(t) * p1**n + k_3 * p1**n)
            #dp2dt = odval * (a_0 * k_2 * profile(t)) / (1 + k_2 * profile(t) + k_3 * p1**n)
            return dp2dt

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

    def residuals(self, df, oddf, rec_df, gamma, plotting=False): 
        def func(x): 
            self.data = []
            self.model = []
            self.concs = []
            a = np.exp(x[0])
            b = np.exp(x[1])
            K = np.exp(x[2])
            n = np.exp(x[3])
            heights = np.exp(x[4:])
            print(np.exp(x))

            p0_1 = 0
            
            residual_array = np.array((0,))
            
            df_sorted = df.sort_values(['Sample', 'Time'])
            oddf_sorted = oddf.sort_values(['Sample', 'Time'])
            #df_sorted = df_sorted[oddf_sorted.Measurement>0]
            #oddf_sorted = oddf_sorted[oddf_sorted.Measurement>0]

            t = df_sorted.Time.unique()
            profile = self.make_profile(heights, t)
            self.profile = profile

            for conc,conc_data in df_sorted.groupby('Concentration1'):
                A = conc
                if np.isnan(A):
                    A = 0
                    odmean = oddf_sorted[oddf_sorted.Concentration1.isna()].groupby('Time').mean()
                    rec_profile = rec_df[rec_df.Concentration1.isna()].sort_values('Time').groupby('Time').mean().Measurement.values
                    rec_profile_t = rec_df[rec_df.Concentration1.isna()].sort_values('Time').groupby('Time').mean().index
                else:
                    odmean = oddf_sorted[oddf_sorted.Concentration1==A].groupby('Time').mean()
                    rec_profile = rec_df[rec_df.Concentration1==A].sort_values('Time').groupby('Time').mean().Measurement.values
                    rec_profile_t = rec_df[rec_df.Concentration1==A].sort_values('Time').groupby('Time').mean().index.values
                
                odval = odmean.Measurement.values
                tt = odmean.index.values
                od = interp1d(tt, odval, bounds_error=False, fill_value='extrapolate')
                sigmean = conc_data.groupby('Time').mean()
                p0_2 = sigmean.Measurement.values[0]
                
                rec_profile = interp1d(rec_profile_t[rec_profile>0], rec_profile[rec_profile>0], bounds_error=False, fill_value='extrapolate')
                
                p,tt = self.forward_model(
                                    a, b, K, n,
                                    od,
                                    profile,
                                    rec_profile,
                                    gamma,
                                    p0_1, p0_2,
                                    tt
                        )
                data = sigmean.Measurement.values
                model = p.ravel()
                self.data.append(data)                    
                self.model.append(model)
                self.concs.append(A)

                if plotting:
                    plt.plot(tt, model, '--')
                    plt.plot(tt, data)
                plt.yscale('log')
                residual_array = np.append(residual_array,  np.log(data[1:]) - np.log(model[1:]))
            print(np.sum(~np.isfinite(residual_array)))
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
            gamma,
            tmin,
            tmax,
            initx,
            bounds
            ):
        # Get biomass time series
        biomass_df = flapjack.analysis(type='Background Correct', 
                            study=study,
                            vector=inverter,
                            media=media,
                            strain=strain,
                            signal=biomass_signal,
                            biomass_signal=biomass_signal,
                            remove_data=False,
                            bg_correction=2
                         )
        biomass_df = biomass_df[(biomass_df.Time>tmin)*(biomass_df.Time<tmax)]

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
                            biomass_signal=biomass_signal,
                            remove_data=False,
                            bg_correction=2)
        rec_df = rec_df[(rec_df.Time>tmin)*(rec_df.Time<tmax)]

        inv_df = flapjack.analysis(vector=inverter,
                            study=study,
                            signal=signal,
                            media=media,
                            strain=strain,
                            type='Background Correct',
                            biomass_signal=biomass_signal,
                            remove_data=False,
                            bg_correction=2)
        inv_df = inv_df[(inv_df.Time>tmin)*(inv_df.Time<tmax)]

        func = self.residuals(inv_df,
                        biomass_df,
                        rec_df,
                        gamma=gamma,
                        plotting=False)
        # Solve for parameters
        res = least_squares(fun=func, 
                            x0=initx,
                            bounds=bounds,
                            diff_step=[1e-1]*len(initx)
                            #ftol=1e-3
        )

        func = self.residuals(inv_df,
                        biomass_df,
                        rec_df,
                        gamma=gamma,
                        plotting=True)
        func(res.x)

        plt.show()

        self.res = res

        print(res)
        xx = np.exp(res.x)
        print(xx)
        self.alpha = xx[:2]
        self.K = xx[2]
        self.n = xx[3]
        heights = xx[4:]
        t = inv_df.Time.unique()
        profile = self.make_profile(heights, t)
        max_profile = profile(t).max()
        self.alpha *= max_profile
        self.profile = self.make_profile(heights, t, normalize=True)
      
