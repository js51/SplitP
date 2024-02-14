"""
Declare constants for the splitp package.
"""

# Currently, this order is needed for the subflattening function 
# to work as it relies on the order of the state space
DNA_state_space = ("A", "C", "G", "T")
DNA_state_space_dict = {state: index for index, state in enumerate(DNA_state_space)}
binary_state_space = (0, 1)