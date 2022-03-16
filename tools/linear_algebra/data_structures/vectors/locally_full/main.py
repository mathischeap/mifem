


from screws.freeze.main import FrozenOnly
from root.config.main import rAnk, mAster_rank, cOmm, np
from tools.linear_algebra.data_structures.vectors.locally_full.do import LocallyFullVectorDo

class LocallyFullVector(FrozenOnly):
    """A LocallyFullVector is a vector: In each core, a full vector is saved. """

    def __init__(self, V):
        """

        :param V:
            - A tuple or list of CSCG forms : we make a LocallyFullVector from their `cochain.globe`
            - GlobalVector or DistributedVector
            - csc_matrix of shape (x, 1)
            - 1d array.
            - int: we make a zero 1d array of shape (V, )

        """
        # ------- parse input ---------------------------------------------------------------
        if isinstance(V, (tuple, list)): # tuple of forms
            if all(hasattr(Vi, 'standard_properties') for Vi in V) and \
                all('CSCG_form' in Vi.standard_properties.tags for Vi in V):

                globe_cochains = list()
                for f in V:
                    if f.cochain.local is not None:
                        globe_cochains.append(LocallyFullVector(f.cochain.globe))
                    else:
                        globe_cochains.append(LocallyFullVector(f.numbering.gathering.GLOBAL_num_dofs))

                v = np.concatenate([f.V for f in globe_cochains])
                self._V_ = v

        elif str(V.__class__.__name__) in ("GlobalVector", "DistributedVector"):
            v = V.V
            v = cOmm.gather(v, root=mAster_rank)
            if rAnk == mAster_rank:
                v = np.sum(v)
            v = cOmm.bcast(v, root=mAster_rank)
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
        assert self._V_.__class__.__name__ == 'ndarray' and np.ndim(self._V_) == 1, "V must be a 1-d array."
        self._do_ = None
        self._freeze_self_()

    @property
    def V(self):
        return self._V_

    @property
    def shape(self):
        return self.V.shape

    def __len__(self):
        return self.shape[0]

    def ___PRIVATE_be_distributed_to___(self, *args, method='sequence'):
        """
        Consider this vector represents cochains of some forms in sequence, we can distribute this vector to the forms.

        :param args: the forms.
        :param method: It can be one of:

            1. ``sequence`` -- We distribute the values in sequence to *args.

        :return:
        """
        if method == 'sequence':
            indices = 0
            V = self.V # V now is 1-d array already.
            for form in args:
                GLOBAL_num_dofs = form.GLOBAL_num_dofs
                # Below we make new locally full vector then distribute it.
                form.cochain.globe = LocallyFullVector(V[indices:indices+GLOBAL_num_dofs])
                indices += GLOBAL_num_dofs
        else:
            raise NotImplementedError(f'distribution method: {method} not coded.')

    @property
    def do(self):
        if self._do_ is None:
            self._do_ = LocallyFullVectorDo(self)
        return self._do_
