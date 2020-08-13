class Metabolism:
    def __init__(self, gamma, mu, dt):
        '''
        Actuator of Cell dynamics
        '''
        self.gamma = gamma
        self.mu = mu
        self.dt = dt
    
    def step(self, proteins):
        self.proteins_np = np.array(proteins)
        step_out = (- self.gamma*self.proteins_np - self.mu*self.proteins_np) * self.dt
        return step_out