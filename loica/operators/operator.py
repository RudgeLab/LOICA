class Operator:
    """
    A class that represents a DNA fragment that encode a genetic operator.

    ...
    
    Attributes
    ----------
    output : Regulator | Reporter
        The output of the operator that is regulated by the input
    uri : str, optional
        SynBioHub URI
    sbol_comp : SBOL Component, optional
        SBOL Component
    name : str, optional
        Name of the operator displayed on the network representation
    color: str, optional
        Color displayed on the network representation
    unit: str, optional
        Units of the characterization data
    """

    def __init__(self, output, name=None, uri=None, sbol_comp=None, color='skyblue', unit='mefl'):
        self.output = output
        self.name = name
        self.uri = uri
        self.sbol_comp = sbol_comp
        self.color = color
        self.unit = unit

    def __str__(self):
        if self.name == None:
            return 'OP'
        else: return self.name