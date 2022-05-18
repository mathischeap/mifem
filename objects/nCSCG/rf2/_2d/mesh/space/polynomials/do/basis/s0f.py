# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/16 9:09 PM
"""
import sys
import numpy as np

if './' not in sys.path: sys.path.append('./')

from objects.nCSCG.rf2._2d.mesh.space.polynomials.do.basis.base import _2nCSCG_RF2_MeshSpacePolynomialsBasisBase


class _2nCSCG_RF2_MeshSpacePolynomialsBasis_S0F(_2nCSCG_RF2_MeshSpacePolynomialsBasisBase):
    """"""

    def __init__(self, mesh, xi_eta):
        """"""
        super(_2nCSCG_RF2_MeshSpacePolynomialsBasis_S0F, self).__init__(mesh, xi_eta)

        if self._xi_eta_.__class__.__name__ == 'Homogeneous':
            assert xi_eta.ndim == 1, f"evaluate_basis only accepts 1-d inputs."
            self._getitem_ = self.___Pr_getitem_Homogeneous___
        elif self._xi_eta_.__class__.__name__ == 'Gauss':
            self._getitem_ = self.___Pr_getitem_Gauss___
        else:
            raise NotImplementedError()
        self._cache_ = dict()
        self._freeze_self_()

    def __getitem__(self, indices):
        """Get the basis functions for the root-cell mesh(indices)"""
        return self._getitem_(indices)

    def ___Pr_getitem_Homogeneous___(self, indices):
        """"""
        level = len(indices)
        cell = self._mesh_(indices)
        N = cell.space.N

        if level == 1:
            if N in self._cache_:
                pass
            else:
                space = cell.space.body
                xi, et = self._xi_eta_[indices]

                bf_xi = space.basises[0].node_basis(x=xi)
                bf_et = space.basises[1].node_basis(x=et)
                bf = np.kron(bf_et, bf_xi)
                _basis_ = (bf,)

                _2d_shape = (len(xi), len(et))

                xi, et = np.meshgrid(xi, et, indexing='ij')
                xi = xi.ravel('F')
                et = et.ravel('F')

                self._cache_[N] = (xi, et), _basis_, _2d_shape

            return self._cache_[N]

        else:
            ind = indices[1:]
            N_ind = (N,) + ind
            if N_ind in self._cache_:
                pass
            else:
                space = cell.space.body
                xi, et = self._xi_eta_[indices]

                bf_xi = space.basises[0].node_basis(x=xi)
                bf_et = space.basises[1].node_basis(x=et)

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

    def ___Pr_getitem_Gauss___(self, indices):
        cell = self._mesh_(indices)
        N = cell.space.N
        if N in self._cache_:
            pass
        else:
            space = cell.space.body
            xi = et = self._xi_eta_[indices][0][0]

            bf_xi = space.basises[0].node_basis(x=xi)
            bf_et = space.basises[1].node_basis(x=et)
            bf = np.kron(bf_et, bf_xi)
            _basis_ = (bf,)

            _2d_shape = (len(xi), len(et))

            xi, et = np.meshgrid(xi, et, indexing='ij')
            xi = xi.ravel('F')
            et = et.ravel('F')

            self._cache_[N] = (xi, et), _basis_, _2d_shape

        return self._cache_[N]







if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
