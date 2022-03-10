"""
Jacobian preconditioner.

"""
from root.config.main import *
from scipy import sparse as spspa
from tools.linear_algebra.preconditioners.base import Preconditioner


class JacobianPreconditioner(Preconditioner):
    """"""
    def __init__(self, A):
        """"""
        super(JacobianPreconditioner, self).__init__(A)
        self._freeze_self_()

    @property
    def M(self):
        A = self._A_.M
        diag = A.diagonal()

        if rAnk != mAster_rank:
            DIAG = None
        else:
            DIAG = np.empty((sIze, self._A_.shape[0]))
        cOmm.Gather(diag, DIAG, root=mAster_rank)

        if rAnk == mAster_rank:
            DIAG = np.sum(DIAG, axis=0)
            DIAG = np.reciprocal(DIAG)
        else:
            DIAG = np.empty((self._A_.shape[0],))

        cOmm.Bcast(DIAG, root=mAster_rank)

        M = spspa.dia_matrix((DIAG, 0), shape=self._A_.shape)

        return M


    @property
    def ___applying_method___(self):
        return 'left_multiply'