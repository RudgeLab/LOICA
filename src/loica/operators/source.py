from .operator import *
import numpy as np
from scipy.optimize import least_squares

class Source(Operator):
    """
    A class that represents a DNA fragment that encode a genetic operator.
    The Source Operator is an abstraction of a constitutive promoter that
    produces output.

    ...
    
    Attributes
    ----------
    output : Regulator | Reporter
        The output of the operator that is constitutively expressed
    rate : float
        Output constitutive expression rate in MEFL/second
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

    def __init__(self, output, rate, uri=None, sbol_comp=None, color='blue', name=None):
        super().__init__(output, name, uri, sbol_comp, color)
        self.rate = rate

    def __str__(self):
        return 'SRC'

    def expression_rate(self, t, dt):
        return self.rate

    def forward_model(
        self,
        Dt=0.25,
        sim_steps=10,
        odval=[1]*97,
        rate=1,
        gamma=0,
        p0=0,
        nt=100
    ):
        p1_list,od_list, A_list,t_list = [],[],[],[]
        p1 = p0
        for t in range(nt):
            p1_list.append(p1)
            t_list.append([t * Dt])
            od = odval[t]
            tt = t*Dt
            for tt in range(sim_steps):
                nextp1 = p1 + (odval[t]*rate - gamma*p1) * Dt / sim_steps
                p1 = nextp1

        ap1 = np.array(p1_list).transpose()
        tt = np.array(t_list).transpose()
        t = np.arange(nt) * Dt
        return ap1,tt

    def residuals(self, df, oddf): 
        def func(x): 
            p0 = x[0]
            rate = x[1]
            residual_list = []
            df_sorted = df.sort_values(['Sample', 'Time'])
            oddf_sorted = oddf.sort_values(['Sample', 'Time'])
            for samp_id,samp_data in df_sorted.groupby('Sample'):
                odval = oddf_sorted[oddf_sorted.Sample==samp_id].Measurement.values
                data = samp_data.Measurement.values
                t = samp_data.Time.values
                dt = np.mean(np.diff(t))
                nt = len(t)
                gamma = 0
                p,tt = self.forward_model(
                            Dt=dt,
                            odval=odval,
                            rate=rate,
                            nt=nt,
                            p0=p0,
                            gamma=gamma
                        )
                model = p[1:]
                residual = (data[1:] - model) 
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
                            ).sort_values(['Sample', 'Time'])

        biomass_df = flapjack.analysis(media=media, 
                            strain=strain,
                            vector=vector,
                            signal=biomass_signal,
                            type='Background Correct',
                            biomass_signal=biomass_signal
                            ).sort_values(['Sample', 'Time'])

        # Bounds for fitting
        lower_bounds = [0, 0]
        upper_bounds = [1e8, 1e8]
        bounds = [lower_bounds, upper_bounds]
        '''
            p0 = x[0]
            rate = x[1]
        '''
        self.residuals_func = self.residuals(
                    expression_df, biomass_df
                    )
        res = least_squares(
                self.residuals_func, 
                [0, 1], 
                bounds=bounds
                )
        self.res = res

        self.p0 = res.x[0]
        self.rate = res.x[1]


