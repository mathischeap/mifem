# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/12 11:57 AM
"""
import sys

import numpy as np

if './' not in sys.path: sys.path.append('./')
from screws.freeze.base import FrozenOnly
from objects.nCSCG.rf2._2d.mesh.space.polynomials.do.basis.allocator import \
    _2nCSCG_RF2_MeshSpacePolynomialsBasisAllocator


class _2nCSCG_SpacePolyDo(FrozenOnly):
    """"""

    def __init__(self, space):
        """"""
        self._space_ = space
        self._freeze_self_()


    def add(self, N):
        """add the space of degree N to the pool."""
        pool = self._space_.___pool___
        assert isinstance(N, int) and N > 0, f"N={N} ({N.__class__.__name__}) is wrong, I need a positive integer."
        assert N not in pool, f"N={N} function space already exists."
        new_space = self._space_.___space_class___((self._space_.___node_type___, N), ndim=2)
        self._space_.___pool___[N] = new_space
        GoN = np.meshgrid(*new_space.nodes, indexing='ij')
        self._space_._GoN_[N] = GoN
        # noinspection PyUnresolvedReferences
        self._GoN_ravel_[N] = GoN[0].ravel('F'), GoN[1].ravel('F')

    def evaluate_basis(self, f, xi_eta):
        """"""
        assert xi_eta.___Pr_is_2nCSCG_RF2_mesh_coo___, f"I need a coo distribution object."
        assert self._space_.mesh.signature == f.signature, f"signature dis-match."
        assert self._space_.mesh.signature == xi_eta.signature, f"signature dis-match."
        basis = _2nCSCG_RF2_MeshSpacePolynomialsBasisAllocator(f.__class__.__name__)(self._space_.mesh, xi_eta)
        return basis




if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
