# -*- coding: utf-8 -*-
from tools.elementwiseCache.gathering.regular.chain_matrix.main import Chain_Gathering_Matrix
from tools.elementwiseCache.dataStructures.operators.bmat.sparse_matrix.helpers.DG import \
    ___BMAT_HELPER_DataGenerator___
from tools.elementwiseCache.dataStructures.operators.bmat.sparse_matrix.helpers.KG import \
    ___BMAT_HELPER_KeyGenerator___


def ___bmat_EWC_sparse_matrices___(blocks, do_a_test=False):
    """A function to do bmat of EWC_SparseMatrix.

    :param blocks:
    :return:
    """
    assert blocks.__class__.__name__ in ('list', 'tuple', 'ndarray'), \
        "please put blocks in list, tuple or array."

    _I = len(blocks)
    J = None
    for i, bR in enumerate(blocks):
        assert bR.__class__.__name__ in ('list', 'tuple', 'ndarray'), \
            "please put blocks in list, tuple or array."
        if J is None:
            J = len(bR)
        else:
            assert J == len(bR)

    elements = None

    RGM = [[None for _ in range(J)] for _ in range(_I)]
    CGM = [[None for _ in range(J)] for _ in range(_I)]

    CLASS = None
    for i, Bi in enumerate(blocks):
        for j, Bij in enumerate(Bi):
            if Bij is None:
                RGM[i][j] = None
                CGM[i][j] = None
            else:

                if CLASS is None:
                    CLASS = Bij.__class__

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

    assert CLASS.__name__ == 'EWC_SparseMatrix', "We must find an EWC_SparseMatrix class."
    EWC = CLASS(elements, DG, KG, bmat_shape=(_I, J))

    # Now we have a look at if we can get the gathering matrices for the bmat result.
    rgm = [None for _ in range(_I)]
    for i in range(_I):
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
        for i in range(_I):
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
    else:
        pass

    if do_a_test and 0 in EWC:
        _ = EWC[0]
        # do a test to check if bmat is fine. this is OKAY even if we will apply
        # some customization later, the cache is done before the customization
    else:
        pass

    return EWC
