# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 11:57 PM
"""
import sys
if './' not in sys.path: sys.path.append('./')

from root.config.main import cOmm, rAnk, mAster_rank, np

from screws.freeze.base import FrozenOnly


class mpRfT2_S0F_Error(FrozenOnly):
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

        F = self._f_.TW.func.___Pr_evaluate_func___()

        detJ = mesh.rcMC.Jacobian(coo)

        local_error = list()
        for rp in mesh.rcfc:
            quad_weights = coo[rp][1][1]
            LEIntermediate =  (v[rp][0] - F[0](*xy[rp]))**n
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

    def distribution(self, visualize=True, density=20, **kwargs):
        """"""
        mesh = self._f_.mesh
        coo = mesh.coo_map.uniform(density, ndim=1)
        xy, v = self._f_.reconstruction(coo, ravel=False)

        f = self._f_.TW.func
        coo = mesh.coo_map.uniform(density, ndim=2)
        __, V = f.reconstruction(coo)
        d = (v - V).magnitude
        if visualize:
            d.visualize(xy, plot_type='contourf', **kwargs)
        return d



if __name__ == '__main__':
    # mpiexec -n 4 python objects/mpRfT/_2d/forms/standard/_0/base/error.py

    from objects.mpRfT._2d.forms.standard._0.inner.main import mpRfT2_Si0F
    from __init__ import rfT2

    mesh = rfT2.rm(100, refinement_intensity=0, N_range=(3,3))

    f = mpRfT2_Si0F(mesh)


    from objects.mpRfT._2d.cf.scalar.main import mpRfT2_Scalar

    def p(t, x, y): return np.sin(np.pi * x) * np.cos(np.pi * y) + t
    s = mpRfT2_Scalar(mesh, p)

    f.TW.func = s
    s.current_time = 0

    f.discretize()

    print(f.error.L())
    f.visualize(show_mesh=True)