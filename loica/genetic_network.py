import networkx as nx

class GeneticNetwork(nx.DiGraph):
    def __init__(self, regulators=None, links=None):
        super().__init__()
        if regulators:
            self.add_nodes_from(regulators)
        if links:
            self.add_edges_from(links)

    def add_tu(self, input_regulator, output_regulator, tu):
        self.add_edge(input_regulator, output_regulator, object=tu)

    def add_regulator(self, regulator):
        self.add_node(regulator)

    def step(circuit, growth_rate=1, dt=0.1):
        for r1,r2,data in circuit.edges.data():
            tu = data['object']
            expression_rate = tu.output(r1.concentration, dt)
            r2.step(expression_rate, growth_rate, dt)
        for r in circuit.nodes:
            r.update()



        
