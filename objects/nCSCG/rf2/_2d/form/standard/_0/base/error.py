# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 11:57 PM
"""
import sys
if './' not in sys.path: sys.path.append('./')

from root.config.main import np, cOmm, rAnk, mAster_rank

from screws.freeze.base import FrozenOnly


class _2nCSCG_RF2_S0F_Error(FrozenOnly):
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
        coo = mesh.coordinates.Gauss(degree_plus)
        xy, v = self._f_.reconstruct(coo, ravel=False)

        F = self._f_.TW.func.___Pr_evaluate_func___()

        local_error = list()
        for i in mesh:
            cell = mesh(i)
            r = repr(cell)
            _ = coo[i][0][1]
            quad_weights = coo[i][1][0]
            detJ = cell.coordinate_transformation.Jacobian(*_)
            LEIntermediate =  (v[r][0] - F[0](*xy[r]))**n

            local_error.append(
                np.einsum('ij, ij, i, j -> ', LEIntermediate, detJ, quad_weights, quad_weights,
                          optimize='optimal'))

        core_local = np.sum(local_error)
        core_local = cOmm.gather(core_local, root=mAster_rank)

        if rAnk == mAster_rank:
            globalError = np.sum(core_local) ** (1 / n)
        else:
            globalError = None
        globalError = cOmm.bcast(globalError, root=mAster_rank)

        return globalError


if __name__ == '__main__':
    # mpiexec -n 4 python objects/nCSCG/rfT2/_2d/form/standard/_0/base/error.py
    from objects.nCSCG.rf2._2d.form.standard._0.inner.main import _2nCSCG_RF2_InnerStandard0Form
    from objects.nCSCG.rf2._2d.__tests__.Random.mesh import random_mesh_of_elements_around as rm2
    mesh = rm2(1000)

    f0 = _2nCSCG_RF2_InnerStandard0Form(mesh)

    from objects.nCSCG.rf2._2d.fields.scalar.main import _2nCSCG_RF2_ScalarField

    def p(t, x, y): return np.sin(np.pi * x) * np.cos(np.pi * y) + t
    s = _2nCSCG_RF2_ScalarField(mesh, p)

    f0.TW.func = s
    s.current_time = 0

    f0.discretize()
    f0.visualize()

    e = f0.error.L()
    print(e)
