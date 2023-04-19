from enum import Enum, auto

class Model(Enum):
    def _generate_next_value_(name, start, count, last_values):
        return name
    JC = auto()
    K2ST = auto()

class Method(Enum):
    def _generate_next_value_(name, start, count, last_values):
        return name
    flattening = auto()
    subflattening = auto()
    distance = auto()
    mutual_information = auto()

class DrawFormat(Enum):
    def _generate_next_value_(name, start, count, last_values):
        return name
    ASCII = auto()
    matplotlib = auto()

class FlatFormat(Enum):
    def _generate_next_value_(name, start, count, last_values):
        return name
    sparse = auto()
    reduced = auto()
