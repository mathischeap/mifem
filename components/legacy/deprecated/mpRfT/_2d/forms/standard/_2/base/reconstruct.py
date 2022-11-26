# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/25 1:42 PM
"""
import sys

if './' not in sys.path: sys.path.append('/')

import numpy as np
from components.freeze.base import FrozenOnly


class mpRfT2_S2F_Reconstruct(FrozenOnly):
    """"""

    def __init__(self, f):
        """"""
        self._f_ = f
        self._freeze_self_()

    def __call__(self, coo_map, rp=None, ravel=False, value_only=False):
        """

        Parameters
        ----------
        coo_map
        rp :
            In which root-cells we are going to reconstruct the standard form.
        ravel

        Returns
        -------

        """
        assert coo_map.___Pr_is_mpRfT2_mesh_coo_map___, f"I need a mpRfT2 coo_map object."
        f = self._f_
        mesh = f.mesh

        # ---- parse indices -----------------------------------------------------------
        if rp is None:
            rps = mesh.rcfc
        else:
            raise NotImplementedError(f"rp={rp} is not implemented.")
        #================================================================================

        if rps is mesh.rcfc:
            full = True
        else:
            full = False

        Basis = mesh.space.do.evaluate_basis(self._f_, coo_map)
        inv_Jac = mesh.rcMC.inverse_Jacobian(Basis)

        if value_only:

            value = dict()

            for rp in rps:
                _, basis, _2d_shape = Basis[rp]
                v_i = np.einsum('ij, i -> j', basis[0]*inv_Jac[rp], f.cochain.local[rp], optimize='greedy')

                if ravel:
                    value[rp] = [v_i, ]
                else:
                    value[rp] = [v_i.reshape(_2d_shape, order='F'), ]

            value = mesh.rcWds.scalar(value, 2, coo_map.distribution, full)

            return value

        else:
            xy = dict()
            value = dict()

            for rp in rps:
                xi_et, basis, _2d_shape = Basis[rp]
                xy_i = mesh[rp].coordinate_transformation.mapping(*xi_et)
                v_i = np.einsum('ij, i -> j', basis[0]*inv_Jac[rp], f.cochain.local[rp], optimize='greedy')

                if ravel:
                    xy[rp] = xy_i
                    value[rp] = [v_i, ]
                else:
                    xy[rp] = [xy_i[0].reshape(_2d_shape, order='F'),
                              xy_i[1].reshape(_2d_shape, order='F')]
                    value[rp] = [v_i.reshape(_2d_shape, order='F'), ]

            xy = mesh.rcWds.vector(xy, 2, coo_map.distribution, full)
            value = mesh.rcWds.scalar(value, 2, coo_map.distribution, full)

            return xy, value

if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
