from collections import UserDict

class Alignment(UserDict):
    def __init__(self, data, taxa):
        self.data = data
        self.taxa = taxa
