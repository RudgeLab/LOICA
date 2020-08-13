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