
class Buffer:
    """
    A class that represents a DNA fragment that encode a genetic operator.
    The Buffer Operator is an abstraction of an inducible promoter that
    maps an input into an output using a Hill function.

    ...
    
    Attributes
    ----------
    input : Regulator | Supplement
        The input of the operator that regulates the expression of the output
    output : Regulator | Reporter
        The output of the operator that is regulated by the input
    a : int | float
        Basal expression rate, LOW
    b : int | float
        Regulated expression rate, HIGH
    K : int | float
        Half expression input concentration
    n : int | float
        Hill coefficient, cooperative degree
    uri : str, optional
        SynBioHub URI
    sbol_comp : SBOL Component, optional
        SBOL Component

    Methods
    -------
    characterize(flapjack, receiver, inverter, media, strain, signal, biomass_signal, gamma)
        Parameterize the Operator model that maps Input concentration into Output expression rate
    """
    def __init__(self, input, output, a, b, K, uri=None, sbol_comp=None):
        self.a = a
        self.b = b
        self.K = K
        self.n = n
        self.input = input
        self.output = output
        self.uri = uri
        self.sbol_comp = sbol_comp
        
    def expression_rate(self, t, dt):
        input_inducer = self.input.concentration
        i = (input_inducer/self.K)**self.n
        expression_rate = ( self.a + self.b*i ) / (1 + i)
        return expression_rate


