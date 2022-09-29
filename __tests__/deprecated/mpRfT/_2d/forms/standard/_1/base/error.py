# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/25 10:22 AM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from root.config.main import cOmm, rAnk, mAster_rank, np

from screws.freeze.base import FrozenOnly


class mpRfT2_S1F_Error(FrozenOnly):
    """"""

    def __init__(self, f):
        """"""
        self._f_ = f
        self._freeze_self_()

    def L(self, n=2, degree_plus=2):
        """

        Parameters
        ----------
        n
        degree_plus

        Returns
        -------

        """
        mesh = self._f_.mesh
        coo = mesh.coo_map.Gauss(degree_plus)
        xy, v = self._f_.reconstruction(coo, ravel=False)

        F = self._f_.analytic_expression.___Pr_evaluate_func___()

        detJ = mesh.rcMC.Jacobian(coo)

        local_error = list()
        for rp in mesh.rcfc:
            quad_weights = coo[rp][1][1]
            LEIntermediate = (v[rp][0] - F[0](*xy[rp])) ** n  + (v[rp][1] - F[1](*xy[rp])) ** n
            local_error.append(
                np.einsum('ij, ij, ij -> ', LEIntermediate, detJ[rp], quad_weights,
                          optimize='optimal'))

        core_local = np.sum(local_error)
        core_local = cOmm.gather(core_local, root=mAster_rank)

        if rAnk == mAster_rank:
            globalError = np.sum(core_local) ** (1 / n)
        else:
            globalError = None
        globalError = cOmm.bcast(globalError, root=mAster_rank)

        return globalError

if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
