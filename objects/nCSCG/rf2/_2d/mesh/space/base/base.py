# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 
"""
import sys

if './' not in sys.path: sys.path.append('./')

from objects.nCSCG.rf2.base.mesh.space.base import nCSCG_SpaceBase


class _2nCSCG_SpaceBase(nCSCG_SpaceBase):
    """"""

    def __init__(self, dN):
        """"""
        super(_2nCSCG_SpaceBase, self).__init__(dN)
        self._ndim_ = 2

    def __repr__(self):
        raise NotImplementedError()




if __name__ == '__main__':
    # mpiexec -n 4 python objects/nCSCG/rfT2/_2d/space/base.py
    pass
