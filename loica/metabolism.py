class Metabolism:
    def __init__(self, gamma, mu, Dt):
        '''
        Actuator of Cell dynamics
        '''
        self.gamma = gamma
        self.mu = mu
        self.Dt = Dt
    
    def step(self, proteins):
        self.proteins_np = np.array(proteins)
        stepOut = (- self.gamma*self.proteins_np - self.mu*self.proteins_np) * self.Dt
        return stepOut