



from tools.linearAlgebra.elementwiseCache.operators.bmat.sparse_matrix.main import \
    ___bmat_EWC_sparse_matrices___




def bmat(blocks):
    """ We bmat some vectors.

    :param blocks:
    :return:
    """
    # check blocks is a list or tuple of lists or tuples and shapes match.
    assert blocks.__class__.__name__ in ('list', 'tuple', 'ndarray'), \
        "please put blocks in list, tuple or array."
    J = None
    for i, bR in enumerate(blocks):
        assert bR.__class__.__name__ in ('list', 'tuple', 'ndarray'), \
            "please put blocks in list, tuple or array."
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
        raise NotImplementedError(
            "Currently, I can only bmat blocks of `EWC_SparseMatrix` and None.")