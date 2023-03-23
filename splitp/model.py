import numpy as np
from scipy.linalg import expm
from . import constants

class model:
    def __init__(self, name=None, state_space=None, init_dist=None, rate_matrix=None):
        """A model of evolution."""
        # Set chosen model properties
        self.name = name if name else "Generic model"
        self.state_space = state_space if state_space else None
        self.init_dist = init_dist if init_dist else None
        self.rate_matrix = rate_matrix if rate_matrix else None
    
    def transition_matrix(self, t):
        """Return the transition matrix for the model."""
        return expm(t * self.rate_matrix)


class GTR(model):
    """
    A general time reversible model of evolution.
    """

    def __init__(
        self,
        name="Generic GTR model",
        state_space=None,
        equilibrium_distribution=None,
        additional_parameters=None
    ):
        self.state_space = state_space
        self.init_dist = equilibrium_distribution
        self.name = name
        self.rate_matrix = None
        n = len(state_space)
        
        ### Error checking
        # Check that we have the correct number of parameters
        if len(additional_parameters) != n * (n - 1) / 2 or len(equilibrium_distribution) != n:
            raise ValueError("Incorrect number of parameters for GTR model.")
        # Check that the equilibrium distribution sums to 1
        if sum(equilibrium_distribution) != 1:
            raise ValueError("Equilibrium distribution must sum to 1.")
        # Check that the parameters are positive
        if any([x < 0 for x in equilibrium_distribution + additional_parameters]):
            raise ValueError("All parameters must be positive.")
        
        ### Construct rate matrix
        # Construct matrix from parameters
        Q = np.zeros((n, n))
        Q[np.triu_indices(n, 1)] = additional_parameters
        Q += Q.T
        # Equilibrium parameters
        E = np.tile(equilibrium_distribution, (n, 1))
        E = E - np.diag(equilibrium_distribution)
        # Construict rate matrix
        Q = np.multiply(Q, E)
        Q -=  np.diag(np.sum(Q, axis=1))
        # Normalise rate matrix
        scale = -np.dot(equilibrium_distribution, np.diag(Q))
        self.rate_matrix = Q / scale

    def __str__(self) -> str:
        return self.name

    @classmethod
    def JukesCantor(cls, rate=1, state_space=constants.DNA_state_space, name="Jukes-Cantor model"):
        """ A Jukes-Cantor model of evolution. """
        n = len(state_space)
        return cls(
            name, 
            state_space, 
            [1/n for _ in range(n)], 
            [rate for _ in range(int(n*(n-1)/2))]
        )

    @classmethod
    def Kimura(cls, transversion_rate=1, transition_rate=1, state_space=constants.DNA_state_space, name="Kimura model"):
        """ A Kimura model of evolution. """
        n = len(state_space)
        # Transitions A-G and C-T
        temp = np.ones((n, n)) * transversion_rate
        A = state_space.index('A')
        G = state_space.index('G')
        C = state_space.index('C')
        T = state_space.index('T')
        temp[A, G] = transition_rate
        temp[G, A] = transition_rate
        temp[C, T] = transition_rate
        temp[T, C] = transition_rate
        parameters = list(temp[np.triu_indices(n, 1)])
        return cls(
            name, 
            state_space, 
            [1/n for _ in range(n)], 
            parameters
        )