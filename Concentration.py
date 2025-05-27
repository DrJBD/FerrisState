class Concentration:
    """
    This class represents a concentration of a substance.
    It can be used to manage the concentration of various substances in a biochemical pathway.
    """
    def __init__(self, substance, concentration):
        self.substance = substance
        self.concentration = concentration

    def __repr__(self):
        return f"{self.substance}: {self.concentration} M"
