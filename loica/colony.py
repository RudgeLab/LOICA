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
        mu = np.exp(-x/self.r0)
        vel = 0.5 * (1 -np.exp(-x/self.r0))
        dx = np.mean(np.diff(x))
        geneprods =  self.circuit.reporters + self.circuit.regulators
        nprods = len(geneprods)
        
        def step(t, y):
            y = y.reshape((len(x),nprods))
            dy = np.zeros_like(y)

            for g,prod in enumerate(geneprods):
                prod.concentration = y[:,g]
                prod.expression_rate = 0

            for o,op in enumerate(self.circuit.operators):
                expression_rate = op.expression_rate(t, 0)
                if type(op.output)==list:
                    for oo in op.output:
                        oo.express(expression_rate)        
                else:
                    op.output.express(expression_rate)

            for g,prod in enumerate(geneprods):
                pi = prod.concentration
                # Finite difference spatial derivatives
                dpidx = np.zeros_like(x)
                for i in range(1, len(x)):
                    dpidx[i] = (pi[i] - pi[i-1]) / dx
                dpidx[0] = pi[0] / dx
                
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
                rkymo[-xx+((t*nx)//nt)-1,t,:] = kymo[xx,t,:]
        return rkymo

    def norm_kymo(self, kymo):
        nkymo = np.zeros_like(kymo)
        for c in range(nkymo.shape[2]):
            nkymo[:,:,c] = kymo[:,:,c] / kymo[:,:,c].max()
        return nkymo

    def kymograph(self, nx, t0, tmax):
        geneprods =  self.circuit.reporters + self.circuit.regulators
        nprods = len(geneprods)
        L = (tmax - t0) / 2
        x = np.linspace(0, L, nx, endpoint=True)

        dx = np.diff(x).mean()
        dt = dx * 0.1
        nt = int((tmax-t0) // dt)
        dydt = self.fun(x)
        y = np.zeros((nx,nprods,nt))
        for g,prod in enumerate(geneprods):
            prod.initialize()
            y[:,g,0] = prod.init_concentration

        for t in range(1,nt):
            y[:,:,t] = y[:,:,t-1] + dydt(t*dt, y[:,:,t-1]).reshape((nx,nprods)) * dt

        nreps = len(self.circuit.reporters)
        kymo = np.zeros((nx,nt,nreps))
        for o in range(nreps):
            kymo[:,:,o] = y[:,o,:]
        return kymo

