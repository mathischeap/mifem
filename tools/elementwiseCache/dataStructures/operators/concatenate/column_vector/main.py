# -*- coding: utf-8 -*-
from tools.elementwiseCache.gathering.regular.chain_matrix.main import Chain_Gathering_Matrix

from tools.elementwiseCache.dataStructures.operators.concatenate.column_vector.helpers.DG import ___concatenate_HELPER_DataGenerator___
from tools.elementwiseCache.dataStructures.operators.concatenate.column_vector.helpers.KG import ___concatenate_HELPER_KeyGenerator___

def ___concatenate_EWC_sparse_vectors___(vectors):
    """"""
    assert vectors.__class__.__name__ in ('list', 'tuple') and len(vectors) > 0, \
        "please put vectors in list or tuple."

    elements = None
    CLASS = None
    for i, v in enumerate(vectors):
        assert v.__class__.__name__ == 'EWC_ColumnVector', f"Cannot handle vectors={vectors}."
        assert v.customize._customizations_ == dict(), \
            f"vectors[{i}] is customized, can not used for concatenate."

        if elements is None:
            CLASS = v.__class__
            elements = v.elements
        else:
            assert elements == v.elements
            assert CLASS == v.__class__

    DG = ___concatenate_HELPER_DataGenerator___(vectors)
    KG = ___concatenate_HELPER_KeyGenerator___(vectors)

    EWC = CLASS(elements, DG, KG, con_shape=(len(vectors),))

    CGMs = list()
    for v in vectors:
        CGMs.append(v.gathering_matrix)

    if None not in CGMs:
        EWC.gathering_matrix = Chain_Gathering_Matrix(CGMs)


    if 0 in EWC: _ = EWC[0]
    # do a test to check if concatenate is fine. This is OKAY even if we will apply
    # some customization later, the cache is done before the customization


    return EWC