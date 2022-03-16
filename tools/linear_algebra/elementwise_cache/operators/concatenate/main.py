"""Store other functions for linear algebra."""

import sys
if './' not in sys.path: sys.path.append('./')

from tools.linear_algebra.elementwise_cache.operators.concatenate.column_vector.main import ___concatenate_EWC_sparse_vectors___
from tools.linear_algebra.elementwise_cache.operators.bmat.main import bmat



def concatenate(vectors):
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
        return ___concatenate_EWC_sparse_vectors___(vectors)
    else:
        raise NotImplementedError(
            "Currently, I can only bmat blocks of `EWC_ColumnVector`.")





if __name__ == '__main__':
    # mpiexec -n 6 python tools/linear_algebra/elementwise_cache/operators/concatenate/main.py
    bmat = bmat
