class Repressor(TU):
    def __init__(self, input, output, a, e, n, dt):
        '''
        IN: Input of repressor
        OUT: Output of the repressor
        '''
        super().__init__(input, output, a, e, n)
        self.dt = dt
        
    def step(self):
        '''Return OUT next state
        '''
        step_out = ((self.a + self.e*(self.input**self.n))/(1 + self.input**self.n)) * self.dt
        return step_out