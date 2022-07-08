class Impactor:
    """
    A class that represents RNA or Protein fragment that encode an impactor.

    ...
    
    Attributes
    ----------
    enzyme : Regulator | Reporter
        The GeneProduct that has active site-impactor
    name : str, optional
        Name of the impactor displayed on the network representation
    color: str, optional
        Color displayed on the network representation

    """

    def __init__(self, enzyme, name=None, color='green'):
        self.enzyme = enzyme
        self.name = name
        self.color = color

    def __str__(self):
        if self.name == None:
            return 'IMPACT'
        else: return self.name