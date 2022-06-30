#  Copyright (c) Conquista Synbio 2022.

import abc

import numpy as np

from monomer import Monomer


class BaseStrand(abc.ABC):
    def __init__(self, content='DNA'):
        self.content = content
        self.strand = None

    @property
    def molecular_weight(self):
        return np.sum(Monomer.weight_table[self.strand])

    def __str__(self):
        return Monomer.conv_table[self.strand]['letter']

    def __len__(self) -> int:
        return len(self.strand)


class Strand(BaseStrand):
    def __init__(self, strand, content: str = 'DNA'):
        super().__init__(content)
        if type(strand) == np.ndarray:
            self.strand = strand
        elif type(strand) == str:
            self.strand = np.array([Monomer.convert(i) for i in strand], dtype=np.ubyte)


class Segment(BaseStrand):
    def __init__(self, start, end, parent: Strand, content='DNA'):
        super().__init__(content)
        self.start = start
        self.end = end

        self.parent = parent

    @property
    def strand(self):
        return self.parent.strand[self.start:self.end]
