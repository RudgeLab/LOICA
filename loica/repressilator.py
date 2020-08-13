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
        self.r1.input = p3
        self.r1.output = p1
        next_p1 = self.r1.step()
        self.r2.input = p1 
        self.r2.output = p2
        next_p2 = self.r2.step()
        self.r3.input = p2
        self.r3.output = p3
        next_p3 = self.r3.step()
        step_out = np.array([next_p1, next_p2, next_p3])
        return step_out