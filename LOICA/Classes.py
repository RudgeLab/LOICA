class Assay:
    def __init__(self, Cell, Environment):
        '''
        Experiment runs a set of cells in parallel, can separate populations in conditions. Fit flapjack
        '''
    pass

class Cell:
    def __init__(self, Metabolism, GeneticNetwork):
        '''
        Cell is a memory of a first order markov desicion process (can be more than one state if needed).
        Coordinate component steps.
        '''
        self.Metabolism  = Metabolism
        self.GeneticNetwork = GeneticNetwork
        self.proteins = self.GeneticNetwork.proteins
        self.memory = [] 
    
    def run(self, steps):
        for _ in range(steps):
            step1 = self.GeneticNetwork.step(self.proteins)
            step2 = self.Metabolism.step(self.proteins)
            Nextproteins = self.proteins + step1 + step2
            self.proteins = Nextproteins
            self.memory.append(self.proteins)

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

class TU:
    def __init__(self, IN, OUT, a, e, n):
        '''
        Transcriptional Unit (TU) is de basic part of this system
        '''
        self.IN = IN
        self.OUT = OUT
        self.a = a
        self.e = e
        self.n = n

class Repressor(TU):
    def __init__(self, IN, OUT, a, e, n, Dt):
        '''
        IN: Input of repressor
        OUT: Output of the repressor
        '''
        super().__init__(IN, OUT, a, e, n)
        self.Dt = Dt
        
    def step(self):
        '''Return OUT next state
        '''
        stepOUT = ((self.a + self.e*(self.IN**self.n))/(1 + self.IN**self.n)) * self.Dt
        #OUT = nextOUT
        return stepOUT

class GeneticNetwork:
    def __init__(self, TU):
        '''
        Actuator of genetic network dynamics
        '''
    pass

class Repressilator:
    '''
    Define the set of proteins that uses
    '''
    
    def __init__(self, proteins, r1, r2, r3):
        self.proteins = proteins
        self.r1 = r1
        self.r2 = r2
        self.r3 = r3
    
    def step(self, proteins):
        '''Think how to randomice the order of reactions if necessary 
        '''
        p1 = proteins[0]
        p2 = proteins[1]
        p3 = proteins[2]
        self.r1.IN = p3
        self.r1.OUT = p1
        Nextp1 = self.r1.step()
        self.r2.IN = p1 
        self.r2.OUT = p2
        Nextp2 = self.r2.step()
        self.r3.IN = p2
        self.r3.OUT = p3
        Nextp3 = self.r3.step()
        stepOut = np.array([Nextp1, Nextp2, Nextp3])
        return stepOut