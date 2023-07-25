class Supplement:
    """
    Representation of a chemical

    ...

    Attributes
    ----------
    name : str
        Name of the supplement
    concentration : int | float
        concentration of the supplement in Molar
    pubchemid : str
        PubChemID URI of the supplement
    supplier_id : str
        Supplier ID of the supplement. An URL of the product that you aquire.
        Accepts list of the form [product URL, catalog number, batch]. 
    sbol_comp : str
        SBOL component of the supplement.
    """
    def __init__(self, name, pubchemid=None, supplier_id=None, sbol_comp=None, color='pink'):
        self.name = name
        self.concentration = 0
        self.pubchemid = pubchemid
        self.color = color
        if type(supplier_id) == list:
            self.supplier_id = supplier_id[0]
            self.supplier_catalog_number = supplier_id[1]
            self.supplier_batch_number = supplier_id[2]
        else:
            self.supplier_id = supplier_id
        self.sbol_comp = sbol_comp

    def __str__(self):
        return self.name

