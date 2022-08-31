# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 8/28/2022 1:58 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from objects.CSCG._3d.discrete_fields.base.main import _3dCSCG_DiscreteField

class _3dCSCG_DF_Vector(_3dCSCG_DiscreteField):
    """"""

    def __init__(self, mesh, coordinates, values, name='no-name'):
        """"""
        super(_3dCSCG_DF_Vector, self).__init__(mesh, coordinates, values, name)
        assert self.vdim == self.mesh.ndim, f"vdim must be 3, now it is {self.vdim}."





if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
