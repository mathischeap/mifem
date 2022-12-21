# -*- coding: utf-8 -*-

from scipy import sparse as spspa


class ___concatenate_HELPER_DataGenerator___:
    """"""
    def __init__(self, vectors):
        self.vectors = vectors
        self._I = len(vectors)

    def __call__(self, i):
        """"""
        output = [None for _ in range(self._I)]
        for j, Vj in enumerate(self.vectors):
            output[j] = Vj[i]
        return spspa.vstack(output, format='csc')
