class model:
    def __init__(self, name=None, state_space=None, init_dist=None, rate_matrix=None):
        """A model of evolution."""

        # Set chosen model properties
        self.name = name if name else "Unnamed Model"
        self.state_space = state_space if state_space else None
        self.init_dist = init_dist if init_dist else None
        self.rate_matrix = rate_matrix if rate_matrix else None
