class Consortium(Sample):
    """
    Representation of microbial consortium (written either to substitute Sample or 
    is a child of sample). Consortium contains Strains, which are connected through 
    extracellular space/medium.
    
    Attributes
    ----------
    strain : Strain
        strain that is part of the consortium
    sample : Sample
        sample to which this consortium belongs
    
    """

    """ 
    I need to have something similar to metabolism here
    So there will be either nested list or dictionary containing wires and their 
    extracelllar concentration

    For example, extracellular_conc = [[wire1, conc1], [wire2, conc2]]
    And each extracellular concentration will be set to 0 in the beginning
    extracellular_conc = [[wire1, 0], [wire2, 0]]

    Since supplement can be a wire, I need to change wire concentration to be equal 
    to supplement

    so let's say supplement = lc.Supplement(name="wire1", concentration=5)
    and wire1.name = "wire1"
    since supplement.name == wire1.name
    extracellular_conc = [[wire1, 5], [wire2, 0] at time 0
    """
    def __init__(self, 
            strain=None, 
            sample=None
            ):
        # need to determine that strains can be list
        self.strain = []
        self.extracellular_concentration = []



    # function that checks whether supplement is the same as wire
    # if yes, set the extracellular concentration of the wire 
    def supplement_is_wire(self, supplement, wire):
        if supplement.name == wire.name:
            extracellular_conc[wire]=supplement.concentration
