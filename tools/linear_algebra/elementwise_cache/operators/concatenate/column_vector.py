


from scipy import sparse as spspa
from tools.linear_algebra.gathering.chain_matrix.main import Chain_Gathering_Matrix
from tools.linear_algebra.elementwise_cache.objects.column_vector.main import EWC_ColumnVector





def ___concatenate_EWC_sparse_vectors___(vectors):
    """"""
    assert isinstance(vectors, (list, tuple)) and len(vectors) > 0, "please put (more than 0) vectors in list or tuple."
    elements = None
    for i, v in enumerate(vectors):
        assert v.__class__.__name__ == 'EWC_ColumnVector', f"Cannot handle vectors={vectors}."
        assert v.customize._customizations_ == dict(), \
            f"vectors[{i}] is customized, can not used for concatenate."

        if elements is None:
            elements = v.elements
        else:
            assert elements == v.elements

    DG = ___concatenate_HELPER_DataGenerator___(vectors)
    KG = ___concatenate_HELPER_KeyGenerator___(vectors)

    EWC = EWC_ColumnVector(elements, DG, KG, con_shape=(len(vectors),))

    CGMs = list()
    for v in vectors:
        CGMs.append(v.gathering_matrix)

    if None not in CGMs:
        EWC.gathering_matrix = Chain_Gathering_Matrix(CGMs)

    return EWC




class ___concatenate_HELPER_DataGenerator___:
    """"""
    def __init__(self, vectors):
        self.vectors = vectors
        self.I = len(vectors)

    def __call__(self, i):
        """"""
        output = [None for _ in range(self.I)]
        for j, Vj in enumerate(self.vectors):
            output[j] = Vj[i]
        return spspa.vstack(output, format='csc')

class ___concatenate_HELPER_KeyGenerator___:
    """"""
    def __init__(self, vectors):
        self.vectors = vectors

    def __call__(self, i):
        """"""
        key_out = ''
        for v in self.vectors:
            key_out += v._KG_(i)
        return key_out