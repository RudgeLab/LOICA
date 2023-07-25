import sbol3
import networkx as nx
from .operators.operator import Operator
from .geneproduct import Regulator, Reporter
from .supplement import Supplement
from .operators.hill1 import Hill1
from .operators.hill2 import Hill2
from .operators.receiver import Receiver
from .operators.source import Source
from typing import List #, Dict, Tuple, Optional, Union, Any
import numpy as np

class GeneticNetwork():
    """
    Representation of a genetic netowrk composed by a set of Operators, Regulators and Reporters.

    ...

    Attributes
    ----------
    operators : List[Operator]
        List of Operators that are part of the genetic network
    regulators : List[Regulator]
        List of Regulators that are part of the genetic network
    reporters : List[Reporter]
        List of Reporters that are part of the genetic network
    vector : int
        Flapjack ID of the vector that is associated with the genetic network

    Methods
    -------
    to_graph()
        Builds a graph representation of the genetic netwok
    draw()
        Generates a plot of the graph representation builded by to_graph()
    to_sbol(sbol_doc=None)
        Generates a SBOL3 Document representation of the genetic network on sbol_doc  
    """
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

    def add_operator(self, ops):
        if issubclass(type(ops), Operator):
            self.operators.append(ops)
        elif type(ops)==list:
            for op in ops:
                if issubclass(type(op), Operator):
                    self.operators.append(op)
                else: print('Unsupported Type, it should be an Operator')
        else: print('Unsupported Type, it should be an Operator')

    def add_regulator(self, regs):
        if issubclass(type(regs), Regulator):
            self.regulators.append(regs)
        elif type(regs)==list:
            for reg in regs:
                if issubclass(type(reg), Regulator):
                    self.regulators.append(reg)
                else: print('Unsupported Type, it should be an Regulator')
        else: print('Unsupported Type, it should be an Regulator')

    def add_reporter(self, reps):
        if issubclass(type(reps), Reporter):
            self.reporters.append(reps)
        elif type(reps)==list:
            for rep in reps:
                if issubclass(type(rep), Reporter):
                    self.reporters.append(rep)
                else: print('Unsupported Type, it should be an Reporter')
        else: print('Unsupported Type, it should be an Reporter')

    def substep_stochastic(self, t=0, dt=0.1, growth_rate=1):
        # Propensities
        a = []

        # Compute expression rates
        for op in self.operators:
            expression_rate = op.expression_rate(t, dt)
            output = op.output
            if type(op.output)!=list:
                output = [output]
            for o in output:
                o.express(expression_rate)

        # Compute propensities for production and degradation of gene products
        gene_products = self.regulators + self.reporters
        for gp in gene_products:
            # Production reeaction
            a.append(gp.expression_rate)
            a.append((gp.degradation_rate + growth_rate) * gp.concentration)

        # Make list of propensities into array
        a = np.array(a)
        # Total of propensities
        A = a.sum()
        
        # Time step
        tau = 1/A * np.log(1/np.random.random())
        # Random number to select next reaction
        a_i = np.random.random() * A

        # Find reaction and update gene product levels
        for i,gp in enumerate(gene_products):
            if a_i < np.sum(a[:i*2+1]):
                # Production of geneproduct gp
                gp.concentration += 1
                break
            elif a_i < np.sum(a[:i*2+2]):
                # Degradation of geneproduct gp
                gp.concentration -= 1
                break

        # Reset expression rates for next step
        for gp in gene_products:
            gp.expression_rate = 0

        # Return elapsed time
        return tau
                
    def step_stochastic(self, growth_rate=1, t=0, dt=0.1):
        delta_t = 0
        while delta_t < dt:
            #print(f'Elapsed time: {delta_t}')
            delta_t += self.substep_stochastic(t=t, dt=dt, growth_rate=growth_rate)
        
    def step(self, growth_rate=1, t=0, dt=0.1):
        for op in self.operators:
            expression_rate = op.expression_rate(t, dt)
            if type(op.output)==list:
                for o in op.output:
                    o.express(expression_rate)
            else:
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
                    for i,inp in enumerate(op.input):
                        type_ = 'positive' if op.alpha[i+1]>op.alpha[0] else 'negative'
                        g.add_edge(inp, op, type=type_)
                else:
                    type_ = 'positive' if op.alpha[1]>op.alpha[0] else 'negative'
                    g.add_edge(op.input, op, type=type_)
            if type(op.output)==list:
                for o in op.output:
                    g.add_edge(op, o, type='positive')
            else:
                g.add_edge(op, op.output, type='positive')
        return g

    def to_contracted_graph(self):
        g = nx.DiGraph()
        for op in self.operators:
            if hasattr(op, 'input'):
                inputs = op.input
                if type(inputs)!=list:
                    inputs = [inputs]
                for i,inp in enumerate(inputs):
                    type_ = 'positive' if op.alpha[i+1]>op.alpha[0] else 'negative'
                    for op2 in self.operators:
                        outputs = op2.output
                        if type(outputs)!=list:
                            outputs = [outputs]
                        if inp in outputs:
                            g.add_edge(op2, op, type=type_)
            outputs = op.output
            if type(outputs)!=list:
                outputs = [outputs]
            for o in outputs:
                if type(o)==Reporter:
                    g.add_edge(op, o, type='positive')
        return g

    def draw(
        self,
        node_shape='o',
        node_size=500,
        linewidths=0,
        alpha=0.5,
        arrowsize=10,
        font_size=6,
        font_family='Tahoma',
        font_weight='bold',
        pos=nx.kamada_kawai_layout,
        contracted=False
        ):
        if contracted:
            g = self.to_contracted_graph()
        else:
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
        neg_edges = [e for e in g.edges(data=True) if e[2]['type']=='negative']
        pos_edges = [e for e in g.edges(data=True) if e[2]['type']=='positive']
        nx.draw_networkx_edges(
            g, 
            pos=pos, 
            width=1, 
            node_shape=node_shape, 
            node_size=node_size,
            arrowsize=arrowsize,
            edgelist=neg_edges,
            connectionstyle='arc3, rad = 0.1',
            arrowstyle='|-|, widthA=0.0, angleA=0, widthB=0.35, angleB=0'
            )
        nx.draw_networkx_edges(
            g, 
            pos=pos, 
            width=1, 
            node_shape=node_shape, 
            node_size=node_size,
            arrowsize=arrowsize,
            edgelist=pos_edges,
            connectionstyle='arc3, rad = 0.1',
            )
        nx.draw_networkx_labels(
            g,
            pos=pos,
            font_size=font_size,
            font_family=font_family,
            font_weight=font_weight
            )

    def to_sbol(self, sbol_doc: sbol3.Document = None) -> sbol3.Document:
        """Convert the genetic network to SBOL.
        :param sbol_doc: The SBOL document to add the genetic network to.
        """
        if sbol_doc:
            doc=sbol_doc
        else: 
            print('No SBOL Document provided')
            print('Generating a new SBOL Document')
            doc = sbol3.Document()
        products = set()
        geneticnetwork = sbol3.Component('geneticnetwork', sbol3.SBO_DNA)
        geneticnetwork.roles.append(sbol3.SO_ENGINEERED_REGION)
        loica_set = set()
        for op in self.operators:
            # Operator Component
            operator_comp = op.sbol_comp
            operator_sc = sbol3.SubComponent(operator_comp)
            # GeneProduct Outputs Component 
            output_str = ''
            output_scs = []
            if type(op.output) != list:
                outputs = [op.output]
            else: 
                outputs = op.output
            for op_output in outputs:
                output_comp = op_output.sbol_comp
                output_scs.append(sbol3.SubComponent(output_comp))
                output_str += f'_{op_output.name}'
            # TODO output string for policistronic operators
            # TU Component
            if type(op)==Source:
                input_str= 'c'
                tu = sbol3.Component(f'TU_{input_str}_{op}{output_str}', sbol3.SBO_DNA) #generalize to multi input/output TUs
                tu.roles.append(sbol3.SO_ENGINEERED_REGION)
                tu.features = [operator_sc] 
                for sc in output_scs:
                    tu.features.append(sc)
            elif type(op)==Hill2: # type(op.input)==List:
                input_str = ''
                for inp in op.input:
                    input_str += f'_{inp.name}'
                tu = sbol3.Component(f'TU{input_str}_{op}{output_str}', sbol3.SBO_DNA) #generalize to multi input/output TUs
                tu.features = [operator_sc] 
                for sc in output_scs:
                    tu.features.append(sc)
                for inp in op.input:
                    input_comp = inp.sbol_comp
                    if type(input_comp)==sbol3.Component:
                        input_sc = sbol3.SubComponent(input_comp)
                        tu.features.append(input_sc)
                    else:
                        tu.features.append(input_comp)
            else:
                input_str= f'_{op.input.name}'
                tu = sbol3.Component(f'TU{input_str}_{op}{output_str}', sbol3.SBO_DNA) #generalize to multi input/output TUs
                tu.features = [operator_sc] 
                for sc in output_scs:
                    tu.features.append(sc)
                input_comp = op.input.sbol_comp
                if type(input_comp)==sbol3.Component:
                    input_sc = sbol3.SubComponent(input_comp)
                    tu.features.append(input_sc)
                else:
                    tu.features.append(input_comp)                  

            tu.roles.append(sbol3.SO_ENGINEERED_REGION)
            for i in range(len(tu.features)-1):
                constraint = sbol3.Constraint(sbol3.SBOL_PRECEDES, tu.features[i], tu.features[i + 1])
                tu.constraints = [constraint]              
            #tu.constraints = [sbol3.Constraint(sbol3.SBOL_PRECEDES, operator_sc, output_sc)]
            # generate a sequence for the TU assuming assembly by type IIS REsnf both parts will have the fusion sites.
            # compare last 4 bp with thefirst 4 bp of the next part, given the preceds constraint.
            # if they are the same then delete one of them and concatenate the rest.
            # else error or comment TU sequence can not be generated, provide ways to add it.

            # Output GeneProduct Component
            for op_output, sc in zip(outputs, output_scs): #make list of tuples? for output_participation
                if type(op_output)==Regulator:
                    if op_output.type_ == 'PRO':
                        output_gp_comp = sbol3.Component(f'{op_output.name}_protein', sbol3.SBO_PROTEIN)
                        output_gp_comp.roles.append(sbol3.SO_TRANSCRIPTION_FACTOR)               
                    elif op_output.type_ == 'RNA':
                        output_gp_comp = sbol3.Component(f'{op_output.name}_rna', sbol3.SBO_RNA)
                        output_gp_comp.roles.append(sbol3.SO_TRANSCRIPTION_FACTOR)
                    else: 
                        print('Unsupported output molecule type')
                elif type(op_output)==Reporter: # For now just support fluorescent reporters
                    if op_output.type_ == 'PRO':
                        output_gp_comp = sbol3.Component(f'{op_output.name}_protein', sbol3.SBO_PROTEIN)
                        output_gp_comp.roles.append('http://purl.obolibrary.org/obo/NCIT_C37894')
                    elif op_output.type_ == 'RNA':
                        output_gp_comp = sbol3.Component(f'{op_output.name}_rna', sbol3.SBO_RNA)
                        output_gp_comp.roles.append('http://purl.obolibrary.org/obo/NCIT_C37894')  
                    else: 
                            print('Unsupported output molecule type')
                else:
                    print('Unsupported output Type')
                output_gp_sc = sbol3.SubComponent(output_gp_comp)
                tu.features.append(output_gp_sc)
                if op_output not in products:
                    products.add(op_output) 
                    loica_set.add(output_gp_comp)
                # Genetic Production Interaction pf the output
                output_participation = sbol3.Participation(roles=[sbol3.SBO_TEMPLATE], participant=sc)
                gp_participation = sbol3.Participation(roles=[sbol3.SBO_PRODUCT], participant=output_gp_sc)
                production = sbol3.Interaction(types=[sbol3.SBO_GENETIC_PRODUCTION], participations=[output_participation, gp_participation])
                tu.interactions.append(production)
            # obtain TU subcomponents sequences, specially CDS and flanking parts sequences
            # look for ATG on the CDS and upstream part sequences (in the case of MoClo the ATG is in the fusion sites)
            # look for stop codons on frame with the ATG.
            # add translated the sequence between the ATG and the stop codon as protein sequence.
            #protein.sequence = tu.cds.sequence

            # Input Product Component
            if type(op) == Source:
                inputs=[]
            elif type(op) == Hill2: #type(op.input) != List:
                inputs = op.input
            else: inputs = [op.input]
            #inputs_prod_sc = []
            for op_input in inputs:
                if type(op_input)==Regulator:
                    if op_input.type_ == 'PRO':
                        input_prod_comp = sbol3.Component(f'{op_input.name}_protein', sbol3.SBO_PROTEIN)
                        input_prod_comp.roles.append(sbol3.SO_TRANSCRIPTION_FACTOR)               
                    elif op_input.type_ == 'RNA':
                        input_prod_comp = sbol3.Component(f'{op_input.name}_rna', sbol3.SBO_RNA)
                        input_prod_comp.roles.append(sbol3.SO_TRANSCRIPTION_FACTOR)
                    else: 
                        print('Unsupported input molecule type')
                elif type(op_input)==Supplement:
                    input_prod_comp = sbol3.Component(f'{op_input.name}_chemical', sbol3.SBO_SIMPLE_CHEMICAL)
                    input_prod_comp.roles.append(sbol3.SO_TRANSCRIPTION_FACTOR)     
                else:
                    print('Unsupported input Type')
                # adds two times prod comp on the repressilator but necessary for normal circuits
                if op_input not in products:
                    products.add(op_input) 
                    loica_set.add(input_prod_comp)
                
                input_prod_sc = sbol3.SubComponent(input_prod_comp)
                tu.features.append(input_prod_sc)
                #inputs_prod_sc.append(input_prod_sc)
                #how can I not create 2 times the same component?
                if type(op_input)!=Regulator: # if it is a regulator it is already created
                    loica_set.add(input_prod_comp)
                 # Input Interaction
                if type(op)==Hill1 and op.alpha[0]>op.alpha[1]:
                    input_participation = sbol3.Participation(roles=[sbol3.SBO_INHIBITOR], participant=input_prod_sc)
                    op_participation = sbol3.Participation(roles=[sbol3.SBO_INHIBITED], participant=operator_sc)
                    interaction = sbol3.Interaction(types=[sbol3.SBO_INHIBITION], participations=[input_participation, op_participation])
                    tu.interactions.append(interaction)
                elif type(op)==Hill2 and op.alpha[0]== max(op.alpha):
                    input_participation = sbol3.Participation(roles=[sbol3.SBO_INHIBITOR], participant=input_prod_sc)
                    op_participation = sbol3.Participation(roles=[sbol3.SBO_INHIBITED], participant=operator_sc)
                    interaction = sbol3.Interaction(types=[sbol3.SBO_INHIBITION], participations=[input_participation, op_participation])
                    tu.interactions.append(interaction)
                elif type(op)==Receiver:
                    input_participation = sbol3.Participation(roles=[sbol3.SBO_STIMULATOR], participant=input_prod_sc)
                    op_participation = sbol3.Participation(roles=[sbol3.SBO_STIMULATED], participant=operator_sc)
                    interaction = sbol3.Interaction(types=[sbol3.SBO_STIMULATION], participations=[input_participation, op_participation])
                    tu.interactions.append(interaction)
                elif type(op)==Hill1 and op.alpha[0]<op.alpha[1]:
                    input_participation = sbol3.Participation(roles=[sbol3.SBO_STIMULATOR], participant=input_prod_sc)
                    op_participation = sbol3.Participation(roles=[sbol3.SBO_STIMULATED], participant=operator_sc)
                    interaction = sbol3.Interaction(types=[sbol3.SBO_STIMULATION], participations=[input_participation, op_participation])
                    tu.interactions.append(interaction)
                elif type(op)==Source:
                    pass
                else:
                    print('Unsupported operator Type')         
            # Model
            #model_string = str(op.__dict__)
            op_model = sbol3.Model(f'LOICA{input_str}_{op}{output_str}_model', 
                            source='https://github.com/SynBioUC/LOICA/blob/master/loica/operators',
                            language='http://identifiers.org/EDAM:format_3996',
                            framework='http://identifiers.org/SBO:0000062',)
                            #attachments=[model_string])
            doc.add(op_model)
            tu.models.append(op_model)
            doc.add(tu)
            tu_sc = sbol3.SubComponent(tu)
            geneticnetwork.features.append(tu_sc)
        loica_list = list(loica_set)
        doc.add(loica_list) 
        if len(geneticnetwork.features) > 1:
            for i in range(len(geneticnetwork.features)-1):
                geneticnetwork.constraints = [sbol3.Constraint(sbol3.SBOL_PRECEDES, geneticnetwork.features[i], geneticnetwork.features[i+1])]
        else: pass
        doc.add(geneticnetwork)
        return doc

        
