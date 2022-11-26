# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/16 6:30 PM
"""
import sys

if './' not in sys.path: sys.path.append('/')

import numpy as np
from components.freeze.base import FrozenOnly


class mpRfT2_S0F_Reconstruct(FrozenOnly):
    """"""

    def __init__(self, f):
        """"""
        self._f_ = f
        self._freeze_self_()

    def __call__(self, coo_map, ravel=False, rp=None, value_only=False):
        """

        Parameters
        ----------
        coo_map
        ravel
        rp :
            In which root-cells we are going to reconstruct the standard form.

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

        if value_only:

            value = dict()

            for rp in rps:
                _, basis, _2d_shape = Basis[rp]
                v_i = np.einsum('ij, i -> j', basis[0], f.cochain.local[rp], optimize='greedy')

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
                v_i = np.einsum('ij, i -> j', basis[0], f.cochain.local[rp], optimize='greedy')

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

    @property
    def matrices(self):
        """The reconstruction matrices."""
        raise NotImplementedError()

    def ___Pr_distribution___(self, sampling_density=20, criterion='L2-norm'):
        """root-cell-wise distribution.

        We use a number to represent the value of the s0f in each root-cell.

        Parameters
        ----------
        sampling_density
        criterion

        Returns
        -------

        """
        mesh = self._f_.mesh

        if criterion == 'L2-norm':
            coo = mesh.coo_map.Gauss(2)
            xy, v = self(coo, value_only=False)

            L2_norm = dict()
            for bc in mesh.basic_cells:
                local_error = list()
                basic_cell = mesh.basic_cells[bc]
                for i in basic_cell:
                    cell = mesh[i]
                    r = repr(cell)
                    _ = coo[i][0][1]
                    quad_weights = coo[i][1][0]
                    detJ = cell.coordinate_transformation.Jacobian(*_)
                    LEIntermediate =  (v[r][0])**2

                    local_error.append(
                        np.einsum('ij, ij, i, j -> ', LEIntermediate, detJ, quad_weights, quad_weights,
                                  optimize='optimal'))

                bc_error = np.sqrt(np.sum(local_error))
                L2_norm[bc] = bc_error

            return L2_norm # dict: keys are local basic-cell number (int), values are the bcW values.

        elif criterion == 'mean-abs':
            coo = mesh.coo_map.uniform(sampling_density, ndim=1)
            v = self(coo, value_only=True)

            mean_abs = dict()
            for bc in mesh.basic_cells:
                local_mean = list()
                basic_cell = mesh.basic_cells[bc]
                for i in basic_cell:
                    cell = mesh[i]
                    r = repr(cell)
                    # noinspection PyTypeChecker
                    local_mean.append(np.mean(np.abs(v[r][0])))

                mean_abs[bc] = np.mean(local_mean)

            return mean_abs # dict: keys are local basic-cell number (int), values are the bcW values.

        else:
            raise NotImplementedError(f"criterion={criterion} not implemented.")







if __name__ == "__main__":
    # mpiexec -n 4 python objects/nCSCG/rfT2/_2d/form/standard/_0/base/reconstruct/main.py
    pass
