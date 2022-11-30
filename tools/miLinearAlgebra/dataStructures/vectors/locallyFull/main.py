# -*- coding: utf-8 -*-

from components.freeze.main import FrozenOnly
from root.config.main import RANK, MASTER_RANK, COMM, np
from tools.miLinearAlgebra.dataStructures.vectors.locallyFull.do import LocallyFullVectorDo
from tools.elementwiseCache.dataStructures.operators.concatenate.main import concatenate


class LocallyFullVector(FrozenOnly):
    """A LocallyFullVector is a vector: In each core, a full vector is saved. """

    def __init__(self, V, chain_method='silly'):
        """

        :param V:
            - A tuple or list of forms : we make a LocallyFullVector from their `cochain.globe`
                When GM is None:
                    we chain the `cochain.globe` in the silly way. Otherwise, we chain them using the
                    GM.
            - GlobalVector or DistributedVector
            - csc_matrix of shape (x, 1)
            - 1d array.
            - int: we make a zero 1d array of shape (V,)
        :param chain_method:
            Only play its role when V is a list (tuple) of forms. It indicates which chain_method
            we use for chaining the dofs of these forms.

        """
        # ------- parse input ---------------------------------------------------------------
        if isinstance(V, (tuple, list)): # tuple of forms
            if all(hasattr(Vi, 'standard_properties') for Vi in V) and \
                all('form' in Vi.standard_properties.tags for Vi in V):

                if chain_method == 'silly':
                    globe_cochains = list()
                    for f in V:
                        if f.cochain.local is not None:
                            globe_cochains.append(LocallyFullVector(f.cochain.globe))
                        else:
                            globe_cochains.append(LocallyFullVector(f.numbering.gathering.global_num_dofs))

                    v = np.concatenate([f.V for f in globe_cochains])
                    self._V_ = v

                elif chain_method == 'sequent':
                    EWC = list()
                    for f in V:
                        if f.cochain.local is not None:
                            EWC.append(f.cochain.EWC)
                        else:
                            from tools.elementwiseCache.dataStructures.objects.columnVector.main import EWC_ColumnVector
                            EWC.append(
                                EWC_ColumnVector(f.mesh.elements, f.num.basis)
                            )
                    EWC = concatenate(EWC)
                    EWC.assembler.chain_method = chain_method
                    EWC.gathering_matrix = V
                    V = EWC.assembled               # GlobalVector
                    V = LocallyFullVector(V).V      # make it a LocallyFullVector
                    self._V_ = V                    # use the V of LocallyFullVector as self._V_

                else:
                    raise NotImplementedError(f"initializing locally full vector from "
                                              f"chain_method={chain_method} "
                                              f"not implemented.")

            else:
                raise NotImplementedError(f"Cannot accept a list: {V}.")
        #___________________________________________________________________________________________
        elif str(V.__class__.__name__) in ("GlobalVector", "DistributedVector"):
            v = V.V
            v = COMM.gather(v, root=MASTER_RANK)
            if RANK == MASTER_RANK:
                v = np.sum(v)
            v = COMM.bcast(v, root=MASTER_RANK)
            self._V_ = v.toarray().ravel()

        elif V.__class__.__name__ == 'csc_matrix':
            assert V.shape[1] == 1, f"V must be of shape (x, 1)."
            self._V_ = V.toarray().ravel()

        elif V.__class__.__name__ == 'ndarray':
            assert np.ndim(V) == 1, f"can only be 1d array."
            self._V_ = V

        elif V.__class__.__name__ in ('int', 'int32', 'int64'):
            self._V_ = np.zeros(V)

        else:
            raise Exception(f"Cannot build LocallyFullVector from {V.__class__.__name__}.")

        # ------- check input ---------------------------------------------------------------
        assert self._V_.__class__.__name__ == 'ndarray' and np.ndim(self._V_) == 1, \
            "V must be a 1-d array."
        self._do_ = None
        self._vector_norm_ = None
        self._freeze_self_()

    def __repr__(self):
        return f"LocallyFullVector{self.shape}:{id(self)}"

    @property
    def V(self):
        return self._V_

    @property
    def shape(self):
        return self.V.shape

    def __len__(self):
        return self.shape[0]

    @property
    def do(self):
        if self._do_ is None:
            self._do_ = LocallyFullVectorDo(self)
        return self._do_

    @property
    def vector_norm(self):
        """"""
        if self._vector_norm_ is None:
            self._vector_norm_ = np.sum(self._V_ ** 2) ** 0.5
        return self._vector_norm_

    def __sub__(self, other):
        """
        LocallyFullVector - LocallyFullVector.

        :param other:
        :return:
        """
        assert other.__class__.__name__ == 'LocallyFullVector', \
            f'Can not do LocallyFullVector - {other.__class__.__name__}'
        assert self.shape == other.shape, f"self shape {self.shape} != other.shape {other.shape}"
        LFV = LocallyFullVector(self.V - other.V)
        return LFV

    def __add__(self, other):
        """
        LocallyFullVector + LocallyFullVector.

        :param other:
        :return:
        """
        assert other.__class__.__name__ == 'LocallyFullVector', \
            f'Can not do LocallyFullVector + {other.__class__.__name__}'
        assert self.shape == other.shape, f"self shape {self.shape} != other.shape {other.shape}"
        LFV = LocallyFullVector(self.V + other.V)
        return LFV