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