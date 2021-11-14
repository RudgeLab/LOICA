import numpy as np
from scipy.optimize import least_squares
from scipy.interpolate import interp1d

class Source:
    def __init__(self, output, rate):
        self.rate = rate
        self.output = output

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

    def residuals(self, data, odval, dt, t): 
        def func(x): 
            nt = len(t)
            p0 = x[0]
            rate = x[1]
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
            residual = (data[1:] - model)  # / tt.ravel()[1:] 
            return residual
        return func

    def characterize(self, flapjack, vector, media, strain, signal, biomass_signal):
        expression_df = flapjack.analysis(media=media, 
                            strain=strain,
                            vector=vector,
                            signal=signal,
                            type='Background Correct',
                            biomass_signal=biomass_signal
                            ).sort_values(['Sample', 'Time'])
        t = expression_df.groupby('Time').mean().index.values
        dt = np.diff(t).mean()
        expression = expression_df.groupby('Time').mean().Measurement.values

        biomass_df = flapjack.analysis(media=media, 
                            strain=strain,
                            vector=vector,
                            signal=biomass_signal,
                            type='Background Correct',
                            biomass_signal=biomass_signal
                            ).sort_values(['Sample', 'Time'])
        biomass = biomass_df.groupby('Time').mean().Measurement.values

        nt = len(t)

        # Bounds for fitting
        lower_bounds = [0, 0]
        upper_bounds = [1e8, 1e8]
        bounds = [lower_bounds, upper_bounds]
        '''
            p0 = x[0]
            rate = x[1]
        '''
        data = expression.ravel()
        self.residuals_func = self.residuals(
                    data, biomass, dt=dt, t=t
                    )
        res = least_squares(
                self.residuals_func, 
                [0, 1], 
                bounds=bounds
                )
        self.res = res

        self.p0 = res.x[0]
        self.rate = res.x[1]


