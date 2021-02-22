import pytest
from splitp import *

@pytest.fixture(scope='class')
def get_splits():
    cases = [None for i in range(0, 9)]
    cases[4+1] = list(tree_helper_functions.all_splits(4, trivial=False))
    return cases

def test_correct_balance(get_splits):
    for i, set in enumerate(get_splits):
        if set != None:
            for split in set:
                if i == 4+1:
                    assert tree_helper_functions.get_balance(split, asTuple=True) == (2,2)
                    assert tree_helper_functions.get_balance(split) == '2|2'