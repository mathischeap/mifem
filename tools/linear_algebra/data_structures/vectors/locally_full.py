



from screws.frozen import FrozenOnly
from root.config import rAnk, mAster_rank, cOmm, np


class LocallyFullVector(FrozenOnly):
    """A LocallyFullVector is a vector: In each core, a full vector is saved. """

    def __init__(self, V):
        if str(V.__class__.__name__) in ("GlobalVector", "DistributedVector"):
            v = V.V
            v = cOmm.gather(v, root=mAster_rank)
            if rAnk == mAster_rank:
                v = np.sum(v)
            v = cOmm.bcast(v, root=mAster_rank)
            self._V_ = v.toarray().ravel()
        elif V.__class__.__name__ == 'csc_matrix':
            self._V_ = V.toarray().ravel()
        elif V.__class__.__name__ == 'ndarray':
            self._V_ = V
        else:
            raise Exception(f"Cannot build LocallyFullVector from {V.__class__.__name__}.")

        assert self._V_ .__class__.__name__ == 'ndarray' and np.ndim(self._V_) == 1, "V must be a 1-d array."

        self._freeze_self_()

    @property
    def V(self):
        return self._V_

    @property
    def shape(self):
        return self.V.shape

    def __len__(self):
        return self.shape[0]

    def DO_distribute_to(self, *args, method='sequence'):
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

