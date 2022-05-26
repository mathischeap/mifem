# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/24 5:13 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from objects.mpRfT._2d.mesh.space.helpers.base import mpRfT2_Mesh_Space_Basis
import numpy as np


class mpRfT2_Mesh_Space_So1F_Basis(mpRfT2_Mesh_Space_Basis):
    def __init__(self, mesh, coo_map):
        """"""
        super(mpRfT2_Mesh_Space_So1F_Basis, self).__init__(mesh, coo_map)
        self._freeze_self_()


    def __getitem__(self, rp):
        """Get the basis functions for the root-cell mesh(indices)"""
        assert isinstance(rp, str)
        return self._getitem_(rp)

    def ___Pr_getitem_uniform___(self, rp):
        """"""
        cell = self._mesh_[rp]
        N = cell.N

        if '-' not in rp:
            if N in self._cache_:
                pass
            else:
                space = cell.space
                xi, et = self._cm_[rp]

                lb_xi = space.basises[0].node_basis(x=xi)
                ed_et = space.basises[1].edge_basis(x=et)
                bf_edge_det = np.kron(ed_et, lb_xi)
                ed_xi = space.basises[0].edge_basis(x=xi)
                lb_et = space.basises[1].node_basis(x=et)
                bf_edge_dxi = np.kron(lb_et, ed_xi)
                _basis_ = (bf_edge_det, bf_edge_dxi)

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

                lb_xi = space.basises[0].node_basis(x=xi)
                ed_et = space.basises[1].edge_basis(x=et)
                if any([lb_xi.shape[1] == 0, ed_et.shape[1] == 0]):
                    bf_edge_det = np.zeros([lb_xi.shape[0] * ed_et.shape[0], 0])
                else:
                    bf_edge_det = np.kron(ed_et, lb_xi)

                ed_xi = space.basises[0].edge_basis(x=xi)
                lb_et = space.basises[1].node_basis(x=et)
                if any([ed_xi.shape[1] == 0, lb_et.shape[1] == 0]):
                    bf_edge_dxi = np.zeros([ed_xi.shape[0] * lb_et.shape[0], 0])
                else:
                    bf_edge_dxi = np.kron(lb_et, ed_xi)

                _basis_ = (bf_edge_det, bf_edge_dxi)

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

            lb_xi = space.basises[0].node_basis(x=xi)
            ed_et = space.basises[1].edge_basis(x=et)
            bf_edge_det = np.kron(ed_et, lb_xi)
            ed_xi = space.basises[0].edge_basis(x=xi)
            lb_et = space.basises[1].node_basis(x=et)
            bf_edge_dxi = np.kron(lb_et, ed_xi)
            _basis_ = (bf_edge_det, bf_edge_dxi)

            _2d_shape = (len(xi), len(et))

            xi, et = np.meshgrid(xi, et, indexing='ij')
            xi = xi.ravel('F')
            et = et.ravel('F')

            self._cache_[N] = (xi, et), _basis_, _2d_shape

        return self._cache_[N]




if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
