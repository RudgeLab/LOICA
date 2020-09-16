import networkx as nx

class GeneticNetwork(nx.DiGraph):
    def __init__(self):
        super().__init__()

    def add_tu(self, tu):
        self.add_node(tu)

    def add_link(self, tu1, tu2, regulator):
        self.add_edge(tu1, tu2, object=regulator)

    def step(self, growth_rate=1, dt=0.1):
        for tu in self.nodes:
            inputs = {reg.name:reg.concentration for tu1,tu2,reg in self.in_edges(tu, data='object')}
            expression_rate = tu.output(inputs, dt)
            for tu1,tu2,regulator in self.out_edges(tu, data='object'):
                regulator.express(expression_rate)

        for tu1,tu2,regulator in self.edges.data('object'):
            regulator.step(growth_rate, dt)



        
