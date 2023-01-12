import collections
import numpy as np
import pandas as pd

def fasta_to_dict(path_to_file):
    file = open(path_to_file, "r")
    alignment = collections.OrderedDict()  # Keeps insertion order
    currentSeq = ""
    for line in file:
        if ">" in line:
            currentSeq = line.replace(">", "").replace("\n", "")
        else:
            if currentSeq in alignment:
                alignment[currentSeq] += line.replace("\n", "")
            else:
                alignment[currentSeq] = line.replace("\n", "")
    if __valid_alignment(alignment):
        return alignment
    else:
        raise ValueError("Invalid alignment")

def read_alignment_from_file(path_to_file):
    file = open(path_to_file, 'r')
    alignment = collections.OrderedDict()  # Keeps insertion order
    currentSeq = ""
    for line in file:
        if '>' in line:
            currentSeq = line.replace('>', '').replace('\n', '')
        else:
            if currentSeq in alignment:
                alignment[currentSeq] += line.replace('\n', '')
            else:
                alignment[currentSeq] = line.replace('\n', '')
    if __valid_alignment(alignment):
        return alignment
    else:
        raise ValueError("Invalid alignment")


def __valid_alignment(alignment):
    ordered = isinstance(
        alignment, collections.OrderedDict
    )  # Must be ordered to avoid incorrect patterns
    sameLength = (
        len(set(len(value) for key, value in alignment.items())) == 1
    )  # Must be same length
    return sameLength and ordered


def get_pattern_counts(alignment):
    patterns = {}
    sequences = list(alignment.values())
    sequenceLength = len(sequences[0])
    usable_sequence_length = 0
    for i in range(sequenceLength):
        pattern_at_site = "".join(s[i] for s in sequences).upper()
        if all(
            c in ("A", "C", "G", "T") for c in pattern_at_site
        ):  # Need to handle U (let U = T)
            usable_sequence_length += 1
            if pattern_at_site in patterns:
                patterns[pattern_at_site] += 1
            else:
                patterns[pattern_at_site] = 1
    return patterns, usable_sequence_length


def pattern_counts_to_probs(patterns, seqLen):
    patterns = patterns.copy()
    for key in patterns.keys():
        patterns[key] = patterns[key] / seqLen
    return patterns


def pattern_probs_from_alignment(path_to_file, return_sequence_length=False):
    alignment = read_alignment_from_file(path_to_file)
    counts, sequence_length = get_pattern_counts(alignment)
    probs = pattern_counts_to_probs(counts, sequence_length)
    if return_sequence_length:
        return probs, sequence_length
    else:
        return probs
