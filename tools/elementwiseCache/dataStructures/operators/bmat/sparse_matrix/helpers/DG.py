# -*- coding: utf-8 -*-
from scipy import sparse as spspa

class ___BMAT_HELPER_DataGenerator___:
    """"""
    def __init__(self, blocks):
        self.blocks = blocks
        self.I = len(blocks)
        self.J = len(blocks[0])

    def __call__(self, item):
        """"""
        output = [[None for _ in range(self.J)] for _ in range(self.I)]
        for i, Bi in enumerate(self.blocks):
            for j, Bij in enumerate(Bi):
                if Bij is not None:
                    output[i][j] = Bij[item]
        return spspa.bmat(output, format='csr')
