#  Copyright (c) Conquista Synbio 2022.

import numpy as np


class Monomer(object):
    dt = np.dtype({'names': ['letter', 'byte'], 'formats': ['S1', np.ubyte]})
    conv_table = np.array([(i, n) for i, n in zip('ATCG', range(4))], dtype=dt)
    weight_table = np.array([i for i in range(4)], dtype=np.float)

    def __init__(self, value: str, content: str = 'DNA'):
        self.value = self.convert(value, content)

    @classmethod
    def convert(cls, s: str, content: str = 'DNA'):
        """ Converts string to number """
        conv = {'DNA': cls.conv_table[:50],
                'RNA': cls.conv_table[50:100],
                'AA': cls.conv_table[100:200],
                }[content]
        conv += conv[200:]
        return conv[conv['letter'] == bytes(s, 'utf-8')]['byte']

    @classmethod
    def revert(cls, n: int):
        """ Converts number to string """
        return cls.conv_table[cls.conv_table['byte'] == n]['letter']

    @property
    def molecular_weight(self) -> float:
        return Monomer.weight_table[self.value]

    def complement(self) -> int:
        """ Generate complementary nucleotide """
        if self.value < 100:
            # A -> T/U; C -> G
            if np.count_nonzero(np.array(0, 2, 50, 52) == self.value):
                return self.value + 1
            # T/U -> A; G -> C
            if np.count_nonzero(np.array(1, 3, 51, 53) == self.value):
                return self.value - 1
        elif np.count_nonzero(np.array(200,201,210,211,256)):
            return self.value

    def __str__(self):
        return self.revert(self.value)
