# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 1:17 AM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly


class _2nCSCG_RF2_MeshBoundaries(FrozenOnly):
    """The boundaries of this 2d nCSCG RF2 mesh."""

    def __init__(self, mesh):
        """"""
        self._mesh_ = mesh
        self._freeze_self_()

    @property
    def names(self):
        """All the mesh boundaries even those having no business with this core."""
        return self._mesh_.cscg.boundaries.names

    @property
    def range_of_segments(self):
        """ dict: Keys are all mesh boundary names, values are the local segments on the boundaries.

        So even there is no local segment on the boundary, the boundary name is still a key of the dict.

        And we have to return it in real time as it may change from time to time.
        """
        #TODO: to be continued.
        raise NotImplementedError()




if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
