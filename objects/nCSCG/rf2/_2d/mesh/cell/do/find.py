# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/08 1:50 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly


class _2nCSCG_CellDoFind(FrozenOnly):
    """"""

    def __init__(self, cell):
        """"""
        self._cell_ = cell
        self._freeze_self_()





if __name__ == "__main__":
    # mpiexec -n 6 python objects/nCSCG/rfT2/_2d/mesh/cell/do/find.py
    # from objects.nCSCG.rfT2._2d.__tests__.Random.mesh import random_mesh_of_elements_around as rm2
    #
    # mesh = rm2(10)
    #
    # from root.save import save
    #
    # save(mesh, 'test_mesh')
    #
    # mesh.visualize()

    from root.read.main import read
    mesh = read('test_mesh.mi')
    # mesh.visualize(show_indices=True)

    if 6 in mesh.cscg.elements:
        cell = mesh((6,))

        cells = cell.do.find.cells_at('U')
        print(cells)


