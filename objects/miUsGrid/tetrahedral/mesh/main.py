# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 9/5/2022 12:30 AM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from objects.miUsGrid.base.mesh.main import miUsGrid_MeshBase


class miUsGrid_TetrahedralMesh(miUsGrid_MeshBase):
    """"""

    def __init__(self):
        """"""
        super(miUsGrid_TetrahedralMesh, self).__init__(3, '3d-name')
        self._freeze_self_()



if __name__ == '__main__':
    # mpiexec -n 4 python objects/miUsGrid/tetrahedral/mesh/main.py
    miUsGrid_TetrahedralMesh()
