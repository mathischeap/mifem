# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/26 5:28 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly


# noinspection PyUnusedLocal
def _bUpper(x, y): return x == 0

# noinspection PyUnusedLocal
def _bDown(x, y): return x == 1

# noinspection PyUnusedLocal
def _bLeft(x, y): return y == 0

# noinspection PyUnusedLocal
def _bRight(x, y): return y == 1


class miUsGrid_TriangularMeshAllocator(FrozenOnly):
    """"""

    def __init__(self):
        """"""
        self._freeze_self_()

    @classmethod
    def ___mesh_path___(cls):
        """"""
        base_path = "./objects/miUsGrid/triangular/mesh/samples/"
        DICT = {
            'test0': base_path + 'test_mesh_0.vtu',
        }
        stK = [2, 4, 6, 8, 10, 12, 14, 16, 32, 64 ,128]
        for key in stK:
            DICT['st' + str(key)] = \
                base_path + f'structured_tests/structured_test_mesh_01_K{key}.vtu'

        return DICT

    @classmethod
    def ___mesh_boundaries___(cls):
        """"""
        boundaries_test0 = {'Upper': _bUpper, 'Down': _bDown, 'Left': _bLeft, 'Right': _bRight}

        DICT = {
            'test0': boundaries_test0,
        }
        stK = [2, 4, 6, 8, 10, 12, 14, 16, 32, 64 ,128]
        for key in stK:
            DICT['st' + str(key)] = boundaries_test0

        return DICT




if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
