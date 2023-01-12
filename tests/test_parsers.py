import pytest
from splitp.parsers import newick, fasta
import os

def test_newick_parser(test_cases_trees):
    for case in test_cases_trees:
        if "newick_string" in case:
            newick_string = case["newick_string"]
            if "JSON_tree" in case:
                JSON_tree = case["JSON_tree"]
                assert newick.newick_to_json(newick_string) == JSON_tree
                assert newick.json_to_newick(JSON_tree) == newick_string
            else:
                assert newick.json_to_newick(newick.newick_to_json(newick_string)) == newick_string

def test_fasta_parser(rootdir):
    test_file = os.path.join(rootdir, 'test_files/test_alignment_1.fa')
    pattern_probs = {
        "ATCG": 2/5,
        "GATC": 1/5,
        "CGAT": 1/5,
        "TCGA": 1/5
    }
    assert fasta.pattern_probs_from_alignment(test_file) == pattern_probs
    counts, length = fasta.get_pattern_counts(fasta.fasta_to_dict(test_file))
    assert fasta.pattern_counts_to_probs(counts, length) == pattern_probs