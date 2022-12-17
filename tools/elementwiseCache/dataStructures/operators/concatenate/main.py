# -*- coding: utf-8 -*-
"""Store other functions for linear algebra.
"""

from tools.elementwiseCache.dataStructures.operators.concatenate.column_vector.main import ___concatenate_EWC_sparse_vectors___

def concatenate(vectors, **kwargs):
    """ We concatenate (vstack) some vectors.

    :param vectors:
    :return:
    """
    assert isinstance(vectors, (list, tuple)), "please put vectors in list or tuple."

    whats_in_vectors = list()
    for v in vectors:
        whats_in_vectors.append(v.__class__.__name__)

    if whats_in_vectors.count('EWC_ColumnVector') == len(whats_in_vectors):
        # vectors are all EWC_ColumnVector
        return ___concatenate_EWC_sparse_vectors___(vectors, **kwargs)
    else:
        raise NotImplementedError(
            "Currently, I can only bmat blocks of `EWC_ColumnVector`.")