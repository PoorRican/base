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

    def complement(self) -> np.ndarray:
        # only perform if there is a contiguous strand
        if np.count_nonzero(self.strand is None):
            raise ValueError('Cannot compute complement for non-contiguous strand')

        complements = np.array([i.complement for i in self.strand])
        return np.invert(complements)

    @abc.abstractmethod
    def insert(self, position: int, seq: __class__) -> None:
        """
        Inserts a str `seq` at `position`

        :returns: None
        """
        raise NotImplementedError

    @abc.abstractmethod
    def delete(self, start: int, end: int) -> None:
        """
        Deletes the str from `start` to `end`

        :returns: None
        """
        raise NotImplementedError

    @abc.abstractmethod
    def modify(self, start: int, end: int, seq: __class__) -> None:
        """
        Replaces the str from `start` to `end` with `seq`

        :returns: None
        """
        raise NotImplementedError


class Strand(BaseStrand):
    def __init__(self, strand, content: str = 'DNA'):
        super().__init__(content)
        if type(strand) == np.ndarray:
            self.strand = strand
        elif type(strand) == str:
            self.strand = np.array([Monomer.convert(i) for i in strand], dtype=np.ubyte)

    def insert(self, position: int, seq: BaseStrand):
        _ = np.concatenate((self.strand[:position], seq.strand))
        self.strand = np.concatenate((_, self.strand[position:]))
        return

    def delete(self, start: int, end: int):
        self.strand = np.concatenate(self.strand[:start], self.strand[end:])
        return

    def modify(self, start: int, end: int, seq: BaseStrand):
        if end - start == len(seq):
            self.strand[start:end] = seq.strand
        else:
            self.delete(start, end)
            self.insert(start, seq)


class Segment(BaseStrand):
    def __init__(self, start, end, parent: Strand, content='DNA'):
        super().__init__(content)
        self.start = start
        self.end = end

        self.parent = parent

    @property
    def strand(self):
        return self.parent.strand[self.start:self.end]

    def insert(self, position: int, seq: BaseStrand):
        return self.parent.insert(position, seq)

    def delete(self, start: int, end: int):
        return self.parent.delete(start, end)

    def modify(self, start: int, end: int, seq: BaseStrand):
        return self.parent.modify(start, end, seq)
