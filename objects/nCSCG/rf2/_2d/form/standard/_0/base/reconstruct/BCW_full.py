# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/16 8:17 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
import numpy as np


class _2nCSCG_RF2_S0F_Reconstruct_BCW_full_LCC(FrozenOnly):
    """"""

    def __init__(self, f):
        """"""
        self._f_ = f
        self._freeze_self_()

    def __call__(self, xi_eta, ravel=False, i=None, value_only=False):
        """

        Parameters
        ----------
        xi_eta
        ravel
        i :
            In which root-cells we are going to reconstruct the standard form.

        Returns
        -------

        """
        assert xi_eta.___Pr_is_2nCSCG_RF2_mesh_coo___, f"I need a coo distribution object."
        assert self._f_.signature == xi_eta.signature, f"signature dis-match."
        f = self._f_
        mesh = f.mesh

        #---- parse i -----------------------------------------------------------
        if i is None:
            INDICES = mesh
        else:
            raise NotImplementedError(f"i={i} is not implemented.")

        if INDICES is mesh:
            full = True
        else:
            full = False


        Basis = mesh.space.do.evaluate_basis(self._f_, xi_eta)

        if value_only:

            value = dict()

            for i in INDICES:
                cell = mesh(i)
                _, basis, _2d_shape = Basis[i]
                v_i = np.einsum('ij, i -> j', basis[0], f.cochain.local[i], optimize='greedy')

                if ravel:
                    value[cell.__repr__()] = [v_i,]
                else:
                    value[cell.__repr__()] = [v_i.reshape(_2d_shape, order='F'), ]

            value = mesh.ids('scalar', value, 2, xi_eta.distribution, full)

            return value

        else:
            xy = dict()
            value = dict()

            for i in INDICES:
                cell = mesh(i)
                xi_et, basis, _2d_shape = Basis[i]
                xy_i = cell.coordinate_transformation.mapping(*xi_et)
                v_i = np.einsum('ij, i -> j', basis[0], f.cochain.local[i], optimize='greedy')

                if ravel:
                    xy[cell.__repr__()] = xy_i
                    value[cell.__repr__()] = [v_i,]
                else:
                    xy[cell.__repr__()] = [xy_i[0].reshape(_2d_shape, order='F'),
                                           xy_i[1].reshape(_2d_shape, order='F')]
                    value[cell.__repr__()] = [v_i.reshape(_2d_shape, order='F'), ]


            xy = mesh.ids('vector', xy, 2, xi_eta.distribution, full)
            value = mesh.ids('scalar', value, 2, xi_eta.distribution, full)

            return xy, value









if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
