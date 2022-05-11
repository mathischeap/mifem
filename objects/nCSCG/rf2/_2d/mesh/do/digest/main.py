# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 12:25 AM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly


class _2nCSCG_MeshDigest(FrozenOnly):
    """"""

    def __init__(self, mesh):
        """"""
        self._mesh_ = mesh
        self._freeze_self_()

    def __call__(self, refinement):
        """We make changes to the mesh according to this refinement object."""
        assert refinement.__class__.__name__ == '_2nCSCG_RF2_Refinement'
        assert refinement.mesh is self._mesh_, f"refinement must possesses me as its mesh."

        #TODO: to be continued.

        self._mesh_.do.update() # update and lock self after digestion!
        raise NotImplementedError()





if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
