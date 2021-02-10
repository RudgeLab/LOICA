from scipy.integrate import solve_ivp
import numpy as np

class Colony:
    def __init__(self, 
            circuit=None, 
            r0=1,
            mu0=1
            ):
        self.circuit = circuit
        self.r0 = r0
        self.mu0 = mu0

    def fun(self, x):
        # Compute growth rates on mesh
        mu = np.exp(-x)
        vel = 0.5 * (1 -np.exp(-x/self.r0))
        dx = np.mean(np.diff(x))
        geneprods = self.circuit.regulators + self.circuit.reporters
        nprods = len(geneprods)
        
        def step(t, y):
            y = y.reshape((len(x),nprods))
            dy = np.zeros_like(y)

            for g,prod in enumerate(geneprods):
                prod.concentration = y[:,g]
                prod.expression_rate = 0

            for o,op in enumerate(self.circuit.operators):
                expression_rate = op.expression_rate(t, 0)
                op.output.express(expression_rate)

            for g,prod in enumerate(geneprods):
                pi = prod.concentration
                # Finite difference spatial derivatives
                dpidx = np.zeros_like(x)
                for i in range(1, len(x)-1):
                    dpidx[i] = (pi[i] - pi[i-1]) / dx
                dpidx[-1] = (pi[-1] - pi[-2]) / dx
                dpidx[0] = (pi[1] - pi[0]) / dx

                # Update protein concs
                dpidt =  prod.expression_rate - (prod.degradation_rate + mu*self.mu0) * pi - vel*dpidx
                dy[:,g] = dpidt

            return dy.ravel()
        return step

    def map_kymo(self, kymo):
        rkymo = np.zeros_like(kymo)
        nx,nt,_ = kymo.shape
        for t in range(nt):
            for xx in range(((t*nx)//nt)):
                rkymo[nx+xx-((t*nx)//nt),t,:] = kymo[xx,t,:]
        return rkymo

    def kymograph(self, nx, nt, t0, tmax):
        geneprods = self.circuit.regulators + self.circuit.reporters
        nprods = len(geneprods)
        L = (tmax - t0) / 2
        x = np.linspace(0, L, nx)
        y0 = np.zeros((nx,nprods))
        for g,prod in enumerate(geneprods):
            prod.initialize()
            y0[:,g] = prod.init_concentration
        y0 = y0.ravel()
        res = solve_ivp(self.fun(x), t_span=(t0,tmax), y0=y0, t_eval=np.linspace(t0,tmax,nt), method='LSODA')
        sol = res.y.reshape((nx,nprods,nt))
        kymo = np.zeros((nx,nt,nprods))
        for o in range(nprods):
            kymo[:,:,o] = sol[:,o,:]
        return self.map_kymo(kymo)


