class Cell:
    def __init__(self, metabolism, genetic_network):
        '''
        Cell is a memory of a first order markov desicion process (can be more than one state if needed).
        Coordinate component steps.
        '''
        self.metabolism  = metabolism
        self.genetic_network = genetic_network
        self.proteins = self.genetic_network.proteins
        self.memory = [] 
    
    def run(self, steps):
        for _ in range(steps):
            step1 = self.genetic_network.step(self.proteins)
            step2 = self.metabolism.step(self.proteins)
            next_proteins = self.proteins + step1 + step2
            self.proteins = next_proteins
            self.memory.append(self.proteins)