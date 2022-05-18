# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 
"""
import sys

if './' not in sys.path: sys.path.append('./')

from objects.nCSCG.rf2.base.mesh.space.base import nCSCG_SpaceBase


class _3nCSCG_SpaceBase(nCSCG_SpaceBase):
    """"""

    def __init__(self, N):
        """"""
        super(_3nCSCG_SpaceBase, self).__init__(N)
        self._ndim_ = 3






if __name__ == '__main__':
    # mpiexec -n 4 python objects/nCSCG/rfT2/_2d/space/base.py
    pass
