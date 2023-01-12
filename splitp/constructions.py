import numpy as np
from splitp.enums import FlatFormat
from scipy.sparse import dok_matrix, coo_matrix
import splitp.constants as constants


def flattening(
    split,
    pattern_probabilities,
    flattening_format=FlatFormat.sparse,
    taxa=None,
):
    if flattening_format is FlatFormat.sparse:
        return __sparse_flattening(split, pattern_probabilities, taxa)
    elif flattening_format is FlatFormat.reduced:
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


def __sparse_flattening(split, pattern_probabilities, taxa):
    format = "dok" # Temporary hard-coded choice
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
            row = __index_of("".join([str(pattern[taxa_indexer[s]]) for s in split[0]]))
            col = __index_of("".join([str(pattern[taxa_indexer[s]]) for s in split[1]]))
            flattening[row, col] = r[1]
        return flattening


def subflattening(split, pattern_probabilities, data=None):
    """
    A faster version of signed sum subflattening. Requires a data dictionary and can be supplied with a bundle of
    re-usable information to reduce the number of calls to the multiplications function.
    """
    state_space = constants.DNA_state_space
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
    sp1, sp2 = len(split[0]), len(split[1])
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
            pattern = __reconstruct_pattern(split, row, col)
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


def __reconstruct_pattern(split, row_label, col_label):
    n = len(split[0]) + len(split[1])
    pattern = {}
    for splindex, loc in enumerate(split[0]):
        pattern[int(str(loc), n) if len(str(loc)) == 1 else int(str(loc)[1:])] = row_label[splindex]
    for splindex, loc in enumerate(split[1]):
        pattern[int(str(loc), n) if len(str(loc)) == 1 else int(str(loc)[1:])] = col_label[splindex]
    return "".join(pattern[i] for i in range(n))
