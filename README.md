<img src="https://user-images.githubusercontent.com/27327007/74098210-a760a700-4b69-11ea-8623-28708864d8c7.png" alt="SplitP" width="200"/>

[![](https://img.shields.io/pypi/v/SplitP.svg)](https://pypi.org/project/SplitP/)  ![](https://github.com/js51/SplitP/workflows/build/badge.svg)
[![Coverage Status](https://coveralls.io/repos/github/js51/SplitP/badge.svg?branch=main&service=github)](https://coveralls.io/github/js51/SplitP?branch=main)

Python package which implements split- and rank-based tools for inferring phylogenies, such as _flattenings_ and _subflattenings_.

### Installation

The latest version of SplitP can be installed via the command
`pip install splitp`

### Examples

Import `splitp` and the associated helper functions
```python
import splitp as sp
from splitp import tree_helper_functions as hf
```
Define trees and work with splits
```python
splits = list(hf.all_splits(4))     # [01|23, 02|13, 03|12]
tree = sp.NXTree('((0,1),(2,3));')	
true_splits = tree.true_splits()    # 01|23
```
Let site patterns evolve under any submodel of the general markov model.
```python
JC_subs_matrix = tree.build_JC_matrix(branch_length:=0.05)
tree.reassign_all_transition_matrices(JC_subs_matrix)
pattern_probs = tree.get_pattern_probabilities()
```
```css
>             0         1
      0    AAAA  0.185844
      1    AAAC  0.003262
      ..    ...       ...
      254  TTTG  0.003262
      255  TTTT  0.185844
```
Simulate sequence alignments from pattern distributions
```python
pattern_frequencies = tree.draw_from_multinomial(pattern_probs, 100)
```
```css
>         0    1
    0  AAAA  0.22
    1  AAAC  0.01
    ..  ...   ...
    2  CCGC  0.03
    3  TTTT  0.14
```
Reconstruct trees using split based methods including flattenings:
```python
F1 = tree.flattening('01|23', pattern_frequencies)
F2 = tree.flattening('02|13', pattern_frequencies)
print(tree.split_score(F1) < tree.split_score(F2))    # True
```
Or subflattenings:
```python
SF = tree.signed_sum_subflattening('01|23', pattern_probs)
print(tree.split_score(SF))   # 0.0
```
For more functionality please see the documentation at [splitp.joshuastevenson.me](http://splitp.joshuastevenson.me/splitp.html).

Please see `CONTRIBUTING.md` for information on contributing to this project.
