<img src="https://user-images.githubusercontent.com/27327007/74098210-a760a700-4b69-11ea-8623-28708864d8c7.png" alt="SplitP" width="200"/>

[![](https://img.shields.io/pypi/v/SplitP.svg)](https://pypi.org/project/SplitP/)  ![](https://github.com/js51/SplitP/workflows/build/badge.svg)
[![Coverage Status](https://coveralls.io/repos/github/js51/SplitP/badge.svg?branch=main&service=github)](https://coveralls.io/github/js51/SplitP?branch=main)
[![DOI](https://zenodo.org/badge/216925092.svg)](https://zenodo.org/badge/latestdoi/216925092)

Python package which implements split- and rank-based tools for inferring phylogenies, such as _flattenings_ and _subflattenings_.

### Installation

The latest version of SplitP can be installed via the command
`pip install splitp`

### Example

The following computes 100 split scores from a 10000bp sequence generated from a balanced 10-taxon tree under the Jukes Cantor model.

```python
# Imports
import splitp

# Set parameters
sequence_length = 10_000
num_taxa = 10
branch_length = 0.05
number_of_splits = 100

# Define model
model = splitp.model.GTR.JukesCantor(1/2)

# Generate tree, splits and alignment
tree = splitp.trees.balanced_newick_tree(num_taxa, branch_length)
splits = splitp.all_splits(tree)
alignment = splitp.generate_alignment(tree, model, sequence_length)

# Compute scores
for s in range(number_of_splits):
    split = next(splits)
    flattening = splitp.flattening(split, alignment, splitp.FlatFormat.reduced)
    score = splitp.split_score(flattening)
    print(split, score)
```
For more functionality please see the documentation at [splitp.joshuastevenson.me](http://splitp.joshuastevenson.me/splitp.html).

Please see `CONTRIBUTING.md` for information on contributing to this project.
