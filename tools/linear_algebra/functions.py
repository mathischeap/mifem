"""
Store other functions for linear algebra.
"""
from tools.linear_algebra.elementwise_cache import ___bmat_EWC_sparse_matrices___
from tools.linear_algebra.elementwise_cache import ___concatenate_EWC_sparse_vectors___


def concatenate(vectors):
    """We concatenate (vstack) some vectors."""
    assert isinstance(vectors, (list, tuple)), "please put vectors in list or tuple."

    whats_in_vectors = list()
    for v in vectors:
        whats_in_vectors.append(v.__class__.__name__)

    if whats_in_vectors.count('EWC_ColumnVector') == len(whats_in_vectors):
        # vectors are all EWC_ColumnVector
        return ___concatenate_EWC_sparse_vectors___(vectors)
    else:
        raise NotImplementedError("Currently, I can only bmat blocks of EWC_ColumnVector.")



def bmat(blocks):
    """

    :param blocks:
    :return:
    """

    # check blocks is a list or tuple of lists or tuples and shapes match.
    assert isinstance(blocks, (list, tuple)), "please put blocks in list or tuple."
    J = None
    for i, bR in enumerate(blocks):
        assert isinstance(bR, (list, tuple)), "please put blocks in list or tuple."
        if J is None:
            J = len(bR)
        else:
            assert J == len(bR), "I need a 2-d list of tuple."

    whats_in_blocks = list()
    for i, Bi in enumerate(blocks):
        for j, Bij in enumerate(Bi):
            if Bij is not None:
                whats_in_blocks.append(Bij.__class__.__name__)
    assert len(whats_in_blocks) > 0, "Only have None in blocks?"

    if whats_in_blocks.count('EWC_SparseMatrix') == len(whats_in_blocks):
        # blocks of EWC_SparseMatrix (and None).
        return ___bmat_EWC_sparse_matrices___(blocks)
    else:
        raise NotImplementedError("Currently, I can only bmat blocks of EWC_SparseMatrix and None.")
