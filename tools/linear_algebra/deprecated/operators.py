# -*- coding: utf-8 -*-
from scipy import sparse as _spspa
import numpy as _np
from tools.linear_algebra.data_structures.global_matrix.main import GlobalMatrix as _GlobalMatrix




def concatenate(vectors, format='csc'):
    a = vectors[0]
    if a.__class__.__name__ == 'GlobalVector':
        for b in vectors[1:]:
            assert b.__class__.__name__ == 'GlobalVector'
    elif a.__class__.__name__ == 'DistributedVector':
        for b in vectors[1:]:
            assert b.__class__.__name__ == 'DistributedVector'
    Vpool = list()
    for v in vectors:
        Vpool.append(v.V)
    V = _spspa.vstack(Vpool, format=format)
    return a.__class__(V)


# noinspection PyShadowingBuiltins
def bmat(blocks, format='csc'):
    """

    :param blocks:
    :param format:
    :return:
    """
    assert _np.ndim(blocks) >= 2, "I need 2-d blocks to bmat them."
    shape = _np.shape(blocks)[:2]
    M_blocks = [[None for _ in range(shape[0])] for __ in range(shape[1])]
    for i, bi in enumerate(blocks):
        for j, bij in enumerate(bi):
            if bij is None:
                M_blocks[i][j] = None
            else:
                assert bij.__class__.__name__ == 'GlobalMatrix', \
                    "I can only bmat GlobalMatrix."
                M_blocks[i][j] = bij.M
    M = _spspa.bmat(M_blocks, format=format)
    return _GlobalMatrix(M)
