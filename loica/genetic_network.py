import sbol3
import networkx as nx
from .geneproduct import Regulator, Reporter

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
            if hasattr(op, 'input'):
                if type(op.input)==list:
                    for i in op.input:
                        g.add_edge(i, op)
                else:
                    g.add_edge(op.input, op)
            g.add_edge(op, op.output)
        return g

    def draw(
        self,
        node_shape='o',
        node_size=200,
        linewidths=2,
        alpha=0.75,
        arrowsize=5,
        font_size=6,
        font_family='Tahoma',
        font_weight='bold',
        pos=nx.kamada_kawai_layout
        ):
        g = self.to_graph()
        pos = pos(g)
        nx.draw_networkx_nodes(
            g, 
            pos=pos, 
            node_color=[n.color for n in g.nodes], 
            node_shape=node_shape, 
            node_size=node_size,
            linewidths=linewidths,
            alpha=alpha
            )
        nx.draw_networkx_edges(
            g, 
            pos=pos, 
            width=1, 
            node_shape=node_shape, 
            node_size=node_size,
            arrowsize=arrowsize
            )
        nx.draw_networkx_labels(
            g,
            pos=pos,
            font_size=font_size,
            font_family=font_family,
            font_weight=font_weight
            )

    def to_sbol(self):
        # prototype for not gate interaction
        doc = sbol3.Document()
        geneticnetwork = sbol3.Component('geneticnetwork', sbol3.SBO_DNA)
        geneticnetwork.roles.append(sbol3.SO_ENGINEERED_REGION)
        for op in self.operators:
            operator_comp = op.sbol_comp
            output_comp = op.output.sbol_comp
            input_comp = op.input.sbol_comp
            operator_sc = sbol3.SubComponent(operator_comp)
            output_sc = sbol3.SubComponent(output_comp)
            input_sc = sbol3.SubComponent(input_comp)
            # TU Component
            tu = sbol3.Component(f'TU_{op.input.name}_{op}_{op.output.name}', sbol3.SBO_DNA) #generalize to multi input/output TUs
            tu.roles.append(sbol3.SO_ENGINEERED_REGION)
            tu.features = [operator_sc, output_sc, input_sc]                  
            tu.constraints = [sbol3.Constraint(sbol3.SBOL_PRECEDES, operator_sc, output_sc)]
            # Interaction
            input_participation = sbol3.Participation(roles=[sbol3.SBO_INHIBITOR], participant=input_comp)
            op_participation = sbol3.Participation(roles=[sbol3.SBO_INHIBITED], participant=operator_comp)
            interaction = sbol3.Interaction(types=[sbol3.SBO_INHIBITION], participations=[input_participation, op_participation])
            tu.interactions.append(interaction)
            # Model
            op_model = sbol3.Model(f'LOICA_{op.input.name}_{op}_{op.output.name}_model', 
                            source='https://github.com/SynBioUC/LOICA/blob/master/loica/operators/not_.py',
                            language='http://identifiers.org/EDAM:format_3996', # http://edamontology.org/format_3996
                            framework='https://github.com/SynBioUC/LOICA'
                            )
            doc.add(op_model)
            tu.models.append(op_model)
            #tu.derived_from = []
            doc.add(tu)
            tu_sc = sbol3.SubComponent(tu)
            geneticnetwork.features.append(tu_sc)
        if len(geneticnetwork.features) > 1:
            for i in range(len(geneticnetwork.features)-1):
                geneticnetwork.constraints = [sbol3.Constraint(sbol3.SBOL_PRECEDES, geneticnetwork.features[i], geneticnetwork.features[i+1])]
        else: pass
        doc.add(geneticnetwork)
        return doc


        
