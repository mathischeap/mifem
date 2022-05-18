# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly

class _2nCSCG_MeshDoFind(FrozenOnly):
    """"""
    def __init__(self, mesh):
        """"""
        self._mesh_ = mesh
        self.___oAd___ = ((-1.,-1.), 2.)
        self._cache_ = dict()
        self._freeze_self_()

    @staticmethod
    def local_indices_of_sub_cell(i):
        """Find the indices of a sub-cell w.r.t. a cell itself.
        """
        if i == 0:
            return 0, 0
        elif i == 1:
            return 1, 0
        elif i == 2:
            return 0, 1
        elif i == 3:
            return 1, 1
        else:
            raise Exception()

    def origin_and_delta(self, indices):
        """consider the level-0 cell (cscg mesh-element) as a [-1,1]^2 reference domain.
        """

        if len(indices) == 1:  # a level-0 cell, namely, a mesh-element of the cscg mesh.
            return self.___oAd___
        else:
            I1 = indices[1:]
            if I1 in self._cache_:
                pass
            else:
                delta = 2 * 0.5 ** (len(I1))
                o = self.origin_and_delta(indices[:-1])[0]
                i, j = self.local_indices_of_sub_cell(indices[-1])
                origin = o[0] + i * delta, o[1] + j * delta
                self._cache_[I1] = origin, delta

            return self._cache_[I1]








if __name__ == '__main__':
    # mpiexec -n 4 python objects/nCSCG/rfT2/_2d/mesh/do/find.py
    from objects.nCSCG.rf2._2d.__tests__.Random.mesh import random_mesh_of_elements_around as rm2
    mesh = rm2(100, refinement_intensity=0.5)
    mesh.visualize()