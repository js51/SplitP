import numpy as np
from splitp.enums import FlatFormat
from scipy.sparse import dok_matrix, coo_matrix
import splitp.constants as constants


def flattening(split, pattern_probabilities, flattening_format=FlatFormat.sparse):
    """
    Compute the flattening of a split given a pattern probability dictionary.

    Args:
        split (str or list): The split to compute the flattening of.
        pattern_probabilities (dict): A dictionary of pattern probabilities.
        flattening_format (FlatFormat): The format to return the flattening in.

    Returns:
        The flattening of the split in the specified format.
    """
    if isinstance(split, str):
        split = split.split("|")
    taxa = sorted(set(split[0]) | set(split[1]))
    if flattening_format is FlatFormat.sparse:
        return __sparse_flattening(split, pattern_probabilities, taxa)
    if flattening_format is FlatFormat.reduced:
        return __reduced_flattening(split, pattern_probabilities, taxa)


def __reduced_flattening(split, pattern_probabilities, taxa):
    if isinstance(split, str):
        split = split.split("|")
    flattening_data = {}
    used_cols = set()
    taxa_indexer = {taxon: i for i, taxon in enumerate(taxa)}
    for r in pattern_probabilities.items():
        pattern = r[0]
        row = __index_of("".join([str(pattern[taxa_indexer[s]]) for s in split[0]]))
        col = __index_of("".join([str(pattern[taxa_indexer[s]]) for s in split[1]]))
        used_cols.add(col)
        try:
            flattening_data[row][col] = r[1]
        except KeyError:
            flattening_data[row] = {col: r[1]}
    column_sort_order = {}

    for i, used_col in enumerate(sorted(used_cols)):
        column_sort_order[used_col] = i

    flattening = np.zeros((len(flattening_data), len(used_cols)))
    for i, (row_index, column_data) in enumerate(sorted(flattening_data.items())):
        for col_index, prob in column_data.items():
            flattening[i, column_sort_order[col_index]] = prob
    return flattening


def __sparse_flattening(
    split, pattern_probabilities, taxa, ban_row_patterns=None, ban_col_patterns=None
):
    format = "dok"  # Temporary hard-coded choice
    if format != "dok":
        raise NotImplementedError("Only dok format is currently supported")
    if isinstance(split, str):
        split = split.split("|")
    taxa_indexer = {taxon: i for i, taxon in enumerate(taxa)}
    if format == "coo":
        rows = []
        cols = []
        data = []
        for r in pattern_probabilities.items():
            if r[1] != 0:
                pattern = r[0]
                row = __index_of(
                    "".join([str(pattern[taxa_indexer[s]]) for s in split[0]])
                )
                col = __index_of(
                    "".join([str(pattern[taxa_indexer[s]]) for s in split[1]])
                )
                rows.append(row)
                cols.append(col)
                data.append(r[1])
        return coo_matrix(
            (data, (rows, cols)), shape=(4 ** len(split[0]), 4 ** len(split[1]))
        )
    elif format == "dok":
        flattening = dok_matrix((4 ** len(split[0]), 4 ** len(split[1])))
        for r in pattern_probabilities.items():
            pattern = r[0]
            row_pattern = "".join([str(pattern[taxa_indexer[s]]) for s in split[0]])
            row = __index_of(row_pattern)
            col_pattern = "".join([str(pattern[taxa_indexer[s]]) for s in split[1]])
            col = __index_of(col_pattern)
            if (ban_col_patterns is not None and col_pattern.count(ban_col_patterns) > 1) or (ban_row_patterns is not None and row_pattern.count(ban_row_patterns) > 1):
                    flattening[row, col] = 0
            else:
                flattening[row, col] = r[1]
        return flattening

sparse_flattening_with_banned_patterns = __sparse_flattening

def subflattening(split, pattern_probabilities, data=None):
    """
    A faster version of signed sum subflattening. Requires a data dictionary and can be supplied with a bundle of
    re-usable information to reduce the number of calls to the multiplications function.
    """
    state_space = constants.DNA_state_space
    taxa = sorted(set(split[0]) | set(split[1]))
    taxa_indexer = {taxon: i for i, taxon in enumerate(taxa)}

    if data is None:
        data = {}
    try:
        coeffs = data["coeffs"]
        labels = data["labels"]
    except KeyError:
        data["coeffs"] = coeffs = {}
        data["labels"] = labels = {}

    if isinstance(split, str):
        split = split.split("|")
    sp1, sp2 = map(len, split)
    subflattening = [[0 for _ in range(3 * sp2 + 1)] for _ in range(3 * sp1 + 1)]
    try:
        row_labels = labels[sp1]
    except KeyError:
        row_labels = list(__subflattening_labels_generator(sp1))
        labels[sp1] = row_labels
    try:
        col_labels = labels[sp2]
    except KeyError:
        col_labels = list(__subflattening_labels_generator(sp2))
        labels[sp2] = col_labels
    banned = (
        {("C", "C"), ("G", "G"), ("A", "T")}
        | {(x, "A") for x in state_space}
        | {("T", x) for x in state_space}
    )
    for r, row in enumerate(row_labels):
        for c, col in enumerate(col_labels):
            pattern = __reconstruct_pattern(split, row, col, taxa_indexer)
            signed_sum = 0
            for table_pattern, value in pattern_probabilities.items():
                try:
                    product = coeffs[(pattern, table_pattern)]
                except KeyError:
                    product = 1
                    for t in zip(pattern, table_pattern):
                        if t not in banned:
                            product *= -1
                    coeffs[(pattern, table_pattern)] = product
                signed_sum += product * value
            subflattening[r][c] = signed_sum
    return np.array(subflattening)


def __index_of(string):
    string = reversed(string)
    index = 0
    for o, s in enumerate(string):
        index += (4**o) * constants.DNA_state_space_dict[s]
    return index


def __subflattening_labels_generator(length):
    n = length
    state_space = constants.DNA_state_space
    other_states = state_space[0:-1]
    special_state = state_space[-1]
    templates = (
        (
            "".join("T" for _ in range(i)),
            "".join("T" for _ in range(n - i - 1)),
        )
        for i in range(n)
    )
    for template in templates:
        for c in other_states:
            yield f"{template[0]}{c}{template[1]}"
    yield "".join(special_state for _ in range(n))


def __reconstruct_pattern(split, row_label, col_label, taxa_indexer):
    n = len(taxa_indexer)
    pattern = {}
    for split_half, label in zip(split, (row_label, col_label)):
        for split_index, taxon in enumerate(split_half):
            pattern[taxa_indexer[taxon]] = label[split_index]
    return "".join(pattern[i] for i in range(n))
