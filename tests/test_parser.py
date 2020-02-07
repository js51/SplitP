import pytest
from SplitP import *

@pytest.fixture(scope='class')
def get_newick_strings():
	cases = []
	cases.append([
		"(((A:0.2, B:0.3):0.3,(C:0.5, D:0.3):0.2):0.3, E:0.7)F:1.0;",
		{
			"id" : "F",
			"children": [
				{"id": "A", "branch_length": 0.1},
				{"id": "B", "branch_length": 0.2},
				{
					"id": "E",
					"branch_length": 0.5,
					"children": [
						{"id": "C", "branch_length": 0.3},
						{"id": "D", "branch_length": 0.4}
						]
				}
			]
		}
	])
	
	
def test_trivial_parsimony(get_trees):
	""" Testing that all trivial splits have parsimony of 1 """
	for case in get_newick_strings:
		assert newick_to_json(case[0]) == case[1]