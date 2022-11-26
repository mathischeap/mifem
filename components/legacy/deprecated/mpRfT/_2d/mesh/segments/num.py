# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/31 4:09 PM
"""
import sys

if './' not in sys.path: sys.path.append('/')

from components.freeze.base import FrozenOnly




class mpRfT2_Mesh_Segments_Num(FrozenOnly):
    """"""

    def __init__(self, segments):
        """"""
        self._segments_ = segments
        self._local_ = 0 # will be updated automatically later.
        self._freeze_self_()

    @property
    def local(self):
        """{int} : amount of local segments"""
        return self._local_

    @property
    def GLOBAL(self):
        """{int} : amount of global segments"""
        raise NotImplementedError()





if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
