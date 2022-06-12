# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/25 2:16 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from objects.mpRfT._2d.mesh.space.helpers.base import mpRfT2_Mesh_Space_Basis
import numpy as np


class mpRfT2_Mesh_Space_S2F_Basis(mpRfT2_Mesh_Space_Basis):
    def __init__(self, mesh, coo_map):
        """"""
        super(mpRfT2_Mesh_Space_S2F_Basis, self).__init__(mesh, coo_map)
        self._freeze_self_()

    @property
    def ___Pr_rcMC_key___(self):
        """"""
        return self._cm_.___Pr_rcMC_key___

    def ___Pr_rcMC_nodes___(self, rp):
        """"""
        return self[rp][0]

    @property
    def ___Pr_sgMC_key___(self):
        """Cannot be used for sgMC"""
        raise Exception("Cannot be used for sgMC")

    def ___Pr_sgMC_nodes___(self, rp):
        """Cannot be used for rcMC"""
        raise Exception("Cannot be used for sgMC")

    def __getitem__(self, rp):
        """Get the basis functions for the root-cell mesh(indices)"""
        assert isinstance(rp, str)
        return self._getitem_(rp)

    def ___Pr_getitem_uniform___(self, rp):
        """"""
        cell = self._mesh_[rp]
        N = cell.N

        if '-' not in rp: # a basic cell
            if N in self._cache_:
                pass
            else:
                space = cell.space
                xi, et = self._cm_[rp]

                bf_xi = space.basises[0].edge_basis(x=xi)
                bf_et = space.basises[1].edge_basis(x=et)
                bf = np.kron(bf_et, bf_xi)
                _basis_ = (bf,)

                _2d_shape = (len(xi), len(et))

                xi, et = np.meshgrid(xi, et, indexing='ij')
                xi = xi.ravel('F')
                et = et.ravel('F')

                self._cache_[N] = (xi, et), _basis_, _2d_shape

            return self._cache_[N]

        else:
            N_ind = str(N) + rp.split('-')[1]

            if N_ind in self._cache_:
                pass
            else:
                space = cell.space
                xi, et = self._cm_[rp]

                bf_xi = space.basises[0].edge_basis(x=xi)
                bf_et = space.basises[1].edge_basis(x=et)

                if any([bf_xi.shape[1] == 0, bf_et.shape[1] == 0]):
                    bf = np.zeros([bf_xi.shape[0] * bf_et.shape[0], 0])
                else:
                    bf = np.kron(bf_et, bf_xi)

                _basis_ = (bf,)

                _2d_shape = (len(xi), len(et))
                xi, et = np.meshgrid(xi, et, indexing='ij')
                xi = xi.ravel('F')
                et = et.ravel('F')

                self._cache_[N_ind] = (xi, et), _basis_, _2d_shape

            return self._cache_[N_ind]

    def ___Pr_getitem_Gauss___(self, rp):
        cell = self._mesh_[rp]
        N = cell.N
        if N in self._cache_:
            pass
        else:
            space = cell.space
            xi = et = self._cm_[rp][0][0]

            bf_xi = space.basises[0].edge_basis(x=xi)
            bf_et = space.basises[1].edge_basis(x=et)
            bf = np.kron(bf_et, bf_xi)
            _basis_ = (bf,)

            _2d_shape = (len(xi), len(et))

            xi, et = np.meshgrid(xi, et, indexing='ij')
            xi = xi.ravel('F')
            et = et.ravel('F')

            self._cache_[N] = (xi, et), _basis_, _2d_shape

        return self._cache_[N]

    def ___Pr_getitem_Lobatto___(self, rp):
        return self.___Pr_getitem_Gauss___(rp)


if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
