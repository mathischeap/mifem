# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 11:23 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.main import FrozenClass


class nCSCG_SpaceBase(FrozenClass):
    """"""

    def __init__(self):
        """"""
        self._mesh_ = None
        self._base_space_ = None
        self._ndim_ = None

    @property
    def mesh(self):
        return self._mesh_

    @mesh.setter
    def mesh(self, mesh):
        if self._mesh_ is None:
            if self.ndim == 2:
                assert mesh.__class__.__name__ == '_2nCSCG_RF2_Mesh'
            elif self.ndim == 3:
                assert mesh.__class__.__name__ == '_3nCSCG_RF2_Mesh'
            else:
                raise Exception()
            self._mesh_ = mesh
        else:
            raise Exception(f"nCSCG space cannot change mesh!")

    @property
    def ndim(self):
        return self._ndim_

    @property
    def base(self):
        return self._base_space_



if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
