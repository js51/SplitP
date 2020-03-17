import pytest
from splitp import *

@pytest.fixture(scope='class')
def get_test_cases():
    cases = []
    matrix_1 = [[0.95, 0.05 / 3, 0.05 / 3, 0.05 / 3], [0.05 / 3, 0.95, 0.05 / 3, 0.05 / 3],
                [0.05 / 3, 0.05 / 3, 0.95, 0.05 / 3], [0.05 / 3, 0.05 / 3, 0.05 / 3, 0.95]]
    t1 = NXTree("(((A,B),C),(D,(E,F)));")
    t1.reassign_all_transition_matrices(np.array(matrix_1))
    cases.append({'tree': t1,
                  'pars_tests': {'012|345': 1, '024|135': 3, '0145|23': 2},
                  })
    return cases

def test_parsimony(get_test_cases):
    """ Testing that all trivial splits have parsimony of 1 """
    for case in get_test_cases:
        trivial_splits = list(all_splits(case['tree'].get_num_taxa(), True, True))
        scores = [case['tree'].parsimony_score(s) for s in trivial_splits]
        corr_triv_scores = all([s == 1 for s in scores])
        pars_tests = case['pars_tests']
        corr_other_scores = all(case['tree'].parsimony_score(test[0]) == test[1] for test in pars_tests.items())
        assert corr_other_scores and corr_triv_scores


def test_subflats_equal(get_test_cases):
    for case in get_test_cases:
        splits = list(all_splits(case['tree'].get_num_taxa(), trivial=False))
        pattern_probs = case['tree'].get_pattern_probabilities()
        for sp in splits[::6]:
            flat = case['tree'].flattening(sp, pattern_probs)
            sf1 = case['tree'].subflattening(case['tree'].transformed_flattening(flat))
            sf2 = case['tree'].subflattening_alt(flat)
            sf3 = case['tree'].sparse_subflattening(sp, pattern_probs).toarray()
            sf4 = case['tree'].signed_sum_subflattening(sp, pattern_probs)
            assert np.all(np.isclose(sf4, sf1))
            assert np.all(np.isclose(sf1, sf2)) 
            assert np.all(np.isclose(sf2, sf3))
