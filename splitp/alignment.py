from collections import UserDict


class Alignment(UserDict):
    def __init__(self, data, taxa):
        self.data = data
        self.taxa = tuple(taxa)
        self.sub_alignments = dict()

    def sub_alignment(self, sub_taxa):
        sub_taxa = tuple(sorted(sub_taxa, key=self.taxa.index))
        try: 
            return self.sub_alignments[sub_taxa]
        except KeyError:
            pass
        sub_taxa_indexer = {taxon: self.taxa.index(taxon) for taxon in sub_taxa}
        sub_taxa_indices = sorted(sub_taxa_indexer.values())
        if not (set(sub_taxa) <= set(self.taxa)):
            raise ValueError("Sub taxa must be a  subset of the taxa")
        all_patterns = set(self.data.keys())
        sub_alignment_data = {}
        while len(all_patterns) > 0:
            pattern = all_patterns.pop()
            trimmed_pattern = "".join(pattern[i] for i in sub_taxa_indices)
            try:
                sub_alignment_data[trimmed_pattern] += self.data[pattern]
            except KeyError:
                sub_alignment_data[trimmed_pattern] = self.data[pattern]
        sub_alignment = Alignment(sub_alignment_data, sub_taxa)
        self.sub_alignments[sub_taxa] = sub_alignment
        return sub_alignment
