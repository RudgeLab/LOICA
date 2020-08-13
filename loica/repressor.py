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