import numpy as np 
from splitp.constructions import flattening, subflattening
from splitp.enums import FlatFormat

taxa = (0,1,2,3)

pattern_probs = {
    "ATCG": 2/5,
    "GATC": 1/5,
    "CGAT": 1/5,
    "TCGA": 1/5
}

splits = [({0,1},{2,3}), ({0,2},{1,3}), ({0,3},{1,2})]
true_split = splits[0]
false_split = splits[1]

indices_to_select = [ 3, 7, 11, 12, 13, 14, 15 ]


def test_flattening():
    F = np.array([
        # AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT
        [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ], #AA 
        [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ], #AC
        [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ], #AG
        [ 0, 0, 0, 0, 0, 0,.4, 0, 0, 0, 0, 0, 0, 0, 0, 0 ], #AT
        [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ], #CA
        [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ], #CC
        [ 0, 0, 0,.2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ], #CG
        [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ], #CT
        [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,.2, 0, 0 ], #GA
        [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ], #GC
        [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ], #GG
        [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ], #GT
        [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ], #TA
        [ 0, 0, 0, 0, 0, 0, 0, 0,.2, 0, 0, 0, 0, 0, 0, 0 ], #TC
        [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ], #TG
        [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ], #TT
        # AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT
    ])

    computed_flattening = flattening(true_split, pattern_probs, taxa=taxa).todense()
    np.testing.assert_array_equal(computed_flattening, F)

    F = np.array([
        # AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT
        [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ], #AA 
        [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,.4, 0 ], #AC
        [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ], #AG
        [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ], #AT
        [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,.2, 0, 0, 0, 0 ], #CA
        [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ], #CC
        [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ], #CG
        [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ], #CT
        [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ], #GA
        [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ], #GC
        [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ], #GG
        [ 0,.2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ], #GT
        [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ], #TA
        [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ], #TC
        [ 0, 0, 0, 0,.2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ], #TG
        [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ], #TT
        # AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT
    ])

    computed_flattening = flattening(false_split, pattern_probs, taxa=taxa).todense()
    np.testing.assert_array_equal(computed_flattening, F)


def test_reduced_flattening():
    F = np.array([
        [  0, .4,  0,  0 ],
        [ .2,  0,  0,  0 ],
        [  0,  0,  0, .2 ],
        [  0,  0, .2,  0 ],
    ])

    computed_flattening = flattening(true_split, pattern_probs, flattening_format=FlatFormat.reduced, taxa=taxa)
    np.testing.assert_array_equal(computed_flattening, F)

    F = np.array([
        [  0,  0,  0, .4 ],
        [  0,  0, .2,  0 ],
        [ .2,  0,  0,  0 ],
        [  0, .2,  0,  0 ],
    ])

    computed_flattening = flattening(false_split, pattern_probs, flattening_format=FlatFormat.reduced, taxa=taxa)
    np.testing.assert_array_equal(computed_flattening, F)

def test_subflattening():
    # Semi-manual computation for the subflattening:
    for split in splits:
        computed_flattening = flattening(split, pattern_probs, taxa=taxa).todense()
        S = np.array([
            [1,-1],
            [1, 1]
        ])
        S0 = np.kron(S, S)
        S = np.kron(S0, S0)    
        transformed_flat = S@computed_flattening@S.T
        subflat = transformed_flat[np.ix_(indices_to_select, indices_to_select)]

        # SplitP computation
        computed_subflattening = subflattening(split, pattern_probs)

        # Test with a tolerance of 1e-15, close enough!
        np.testing.assert_allclose(computed_subflattening, subflat, rtol=1e-15, atol=1e-15)