# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 11:23 PM
"""
import sys
if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly



class nCSCG_SpaceBase(FrozenOnly):
    """"""
    def __init__(self, dN):
        """"""
        assert isinstance(dN, int) and dN > 0, f"dN={dN} wrong, must be a positive integer."
        self._PRM = None # parameters; for restoring the space.
        self.___pool___ = dict() # the pool of all useful function spaces. keys are their degrees.
        self.___space_class___ = None # the cscg space class. All spaces will be of this type.

        self._dN_ = dN # default N
        self._ndim_ = None
        self._mesh_ = None



    def __repr__(self):
        raise NotImplementedError()

    def __contains__(self, N):
        return N in self.___pool___

    def __iter__(self):
        for N in self.___pool___:
            yield N

    def __len__(self):
        return len(self.___pool___)

    def __getitem__(self, N):
        return self.___pool___[N]



    @property
    def dN(self):
        """the basis function degree of the default function space."""
        return self._dN_

    @property
    def ndim(self):
        return self._ndim_

    @property
    def mesh(self):
        return self._mesh_










if __name__ == '__main__':
    # mpiexec -n 4 python objects/nCSCG/rfT2/base/space/main.py
    pass
