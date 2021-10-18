import networkx as nx
from .geneproduct import Regulator

class GeneticNetwork():
    def __init__(self, vector=None):
        self.operators = []
        self.regulators = []
        self.reporters = []
        self.vector = vector

    def initialize(self):
        for regulator in self.regulators:
            regulator.initialize()
        for reporter in self.reporters:
            reporter.initialize()

    def add_operator(self, op):
        self.operators.append(op)

    def add_regulator(self, reg):
        self.regulators.append(reg)

    def add_reporter(self, rep):
        self.reporters.append(rep)

    def add_operators(self, ops):
        for op in ops:
            self.operators.append(op)

    def add_regulators(self, regs):
        for reg in regs:
            self.regulators.append(reg)

    def add_reporters(self, reps):
        for rep in reps:
            self.reporters.append(rep)

    def step(self, growth_rate=1, t=0, dt=0.1):
        for op in self.operators:
            expression_rate = op.expression_rate(t, dt)
            op.output.express(expression_rate)

        for regulator in self.regulators:
            regulator.step(growth_rate, dt)

        for reporter in self.reporters:
            reporter.step(growth_rate, dt)

    def to_graph(self):
        g = nx.DiGraph()
        for op in self.operators:
            if type(op.output)==Regulator:
                if type(op.input)==list:
                    for i in op.input:
                        g.add_edge(i, op)
                else:
                    g.add_edge(op.input, op)
                g.add_edge(op, op.output)
        return g

    def draw(self, node_shape='s', node_size=600, pos=nx.circular_layout):
        g = self.to_graph()
        pos = pos(g)
        nx.draw_networkx_nodes(
            g, 
            pos=pos, 
            node_color=[n.color for n in g.nodes], 
            node_shape=node_shape, 
            node_size=node_size
            )
        nx.draw_networkx_edges(
            g, 
            pos=pos, 
            width=1, 
            node_shape=node_shape, 
            node_size=node_size
            )
        nx.draw_networkx_labels(g, pos=nx.circular_layout(g))

    def to_sbol(self):
        #Operators
        operators = []
        for op in self.operators:
            operator_doc = op.sbol_doc
            geneproduct_doc = op.output.sbol_doc
            operator_sc = sbol3.SubComponent(operator_doc)
            geneproduct_sc = sbol3.SubComponent(geneproduct_doc)

            tu = sbol3.Component('tu', sbol3.SBO_DNA)
            tu.roles.append(sbol3.SO_ENGINEERED_REGION)
            tu.features = [operator_sc, geneproduct_sc]
            tu.constraints = [sbol3.Constraint(sbol3.SBOL_PRECEDES, operator_sc, geneproduct_sc)]

            doc = sbol3.Document()
            doc.add(operator_doc)
            doc.add(geneproduct_doc)
            doc.add(tu)
            operators.append(doc)
        return operators


        
