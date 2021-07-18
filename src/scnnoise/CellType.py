class cellType:
    
    def __init__(self, rxn_rates, children):
        """
        Takes in rxn_rates and children and constructs a cellType node
        """
        self.rxn_rates = rxn_rates
        self.children = children
