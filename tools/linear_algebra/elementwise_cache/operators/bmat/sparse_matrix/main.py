


from tools.linear_algebra.gathering.chain_matrix.main import Chain_Gathering_Matrix
from tools.linear_algebra.elementwise_cache.objects.sparse_matrix.main import EWC_SparseMatrix
from tools.linear_algebra.elementwise_cache.operators.bmat.sparse_matrix.helpers.DG import ___BMAT_HELPER_DataGenerator___
from tools.linear_algebra.elementwise_cache.operators.bmat.sparse_matrix.helpers.KG import ___BMAT_HELPER_KeyGenerator___





def ___bmat_EWC_sparse_matrices___(blocks):
    """A function to do bmat of EWC_SparseMatrix.

    :param blocks:
    :return:
    """
    assert isinstance(blocks, (list, tuple)), "please put blocks in list or tuple."

    I = len(blocks)
    J = None
    for i, bR in enumerate(blocks):
        assert isinstance(bR, (list, tuple)), "please put blocks in list or tuple."
        if J is None:
            J = len(bR)
        else:
            assert J == len(bR)

    elements = None

    RGM = [[None for _ in range(J)] for _ in range(I)]
    CGM = [[None for _ in range(J)] for _ in range(I)]

    for i, Bi in enumerate(blocks):
        for j, Bij in enumerate(Bi):
            if Bij is None:
                RGM[i][j] = None
                CGM[i][j] = None
            else:
                assert Bij.customize._customizations_ == dict(), \
                    f"customized block[{i}][{j}] cannot be used for bmat."
                # this is because of the cache function. the customizations of block will not be renewed in cache.

                assert Bij.__class__.__name__ == 'EWC_SparseMatrix', \
                    f"I only handle EWC_SparseMatrix, but block[{i}][{j}] now is {Bij.__class__.__name__}."

                if elements is None:
                    elements = Bij.elements
                else:
                    assert Bij.elements == elements

                RGM[i][j], CGM[i][j] = Bij.gathering_matrices

    assert elements is not None, f"blocks of only None?"
    DG = ___BMAT_HELPER_DataGenerator___(blocks)
    KG = ___BMAT_HELPER_KeyGenerator___(blocks)
    EWC = EWC_SparseMatrix(elements, DG, KG, bmat_shape=(I, J))

    # Now we have a look at if we can get the gathering matrices for the bmat result.
    rgm = [None for _ in range(I)]
    for i in range(I):
        for j in range(J):
            Rij = RGM[i][j]
            if Rij is None:
                pass
            else:
                # noinspection PyUnresolvedReferences
                assert Rij.__class__.__name__ == "Chain_Gathering_Matrix"
                if rgm[i] is None:
                    rgm[i] = Rij
                else:
                    assert rgm[i] == Rij, f"bmat wrong blocks. Gathering matrices for row[{i}] dis-match."

    cgm = [None for _ in range(J)]
    for j in range(J):
        for i in range(I):
            Cij = CGM[i][j]
            if Cij is None:
                pass
            else:
                # noinspection PyUnresolvedReferences
                assert Cij.__class__.__name__ == "Chain_Gathering_Matrix"
                if cgm[j] is None:
                    cgm[j] = Cij
                else:
                    assert cgm[j] == Cij, f"bmat wrong blocks. Gathering matrices for col[{j}] dis-match."

    doR = False if None in rgm else True
    doC = False if None in cgm else True

    if doR and doC:
        R_CGM = Chain_Gathering_Matrix(rgm)
        C_CGM = Chain_Gathering_Matrix(cgm)
        EWC.gathering_matrices = (R_CGM, C_CGM)

    if 0 in EWC: _ = EWC[0]
    # do a test to check if bmat is fine. this is OKAY even if we will apply
    # some customization later, the cache is done before the customization

    return EWC
