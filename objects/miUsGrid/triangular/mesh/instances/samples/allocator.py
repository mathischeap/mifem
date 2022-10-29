# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/26 5:28 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from objects.miUsGrid.triangular.mesh.instances.samples.rand_mesh_0 import generate_rand_mesh_0
from root.config.main import COMM


# noinspection PyUnusedLocal
def _bUpper(x, y): return x == 0

# noinspection PyUnusedLocal
def _bDown(x, y): return x == 1

# noinspection PyUnusedLocal
def _bLeft(x, y): return y == 0

# noinspection PyUnusedLocal
def _bRight(x, y): return y == 1

# noinspection PyUnusedLocal
def __bUpper__(x, y): return x == -1

# noinspection PyUnusedLocal
def __bDown__(x, y): return x == 1

# noinspection PyUnusedLocal
def __bLeft__(x, y): return y == -1

# noinspection PyUnusedLocal
def __bRight__(x, y): return y == 1




class miUsGrid_TriangularMeshAllocator(FrozenOnly):
    """"""

    def __init__(self):
        """"""
        self._freeze_self_()

    @classmethod
    def ___mesh_path___(cls):
        """"""
        COMM.barrier()
        generate_rand_mesh_0()
        base_path = str(__file__).split('allocator')[0]

        #------- mesh examples ---------------------------------------------------------------------
        DICT = {
            'test0': base_path + 'test_mesh_0.vtu',
            'test2': base_path + 'test2.vtu',
            'rand0': base_path + 'rand0.vtu'
        }
        #------- stK meshes ------------------------------------------------------------------------
        stK = [2, 4, 6, 8, 10, 12, 14, 16, 32, 48, 64 ,128]
        for key in stK:
            DICT['st' + str(key)] = \
                base_path + f'structured_tests/structured_test_mesh_01_K{key}.vtu'
        #------- stKd meshes -----------------------------------------------------------------------
        DICT.update(
            {
                'sqK8d2': base_path + f'structured_square/K8d2.vtu',
                'sqK16d2': base_path + f'structured_square/K16d2.vtu',
                'sqK32d2': base_path + f'structured_square/K32d2.vtu',
                'sqK48d2': base_path + f'structured_square/K48d2.vtu',
                'sqK64d2': base_path + f'structured_square/K64d2.vtu',
                'sqK128d2': base_path + f'structured_square/K128d2.vtu',
            }
        )
        COMM.barrier()
        #___________________________________________________________________________________________
        return DICT

    @classmethod
    def ___mesh_boundaries___(cls):
        """"""
        boundaries_test0 = {'Upper': _bUpper, 'Down': _bDown, 'Left': _bLeft, 'Right': _bRight}
        boundaries_test1 = {'Upper': __bUpper__, 'Down': __bDown__, 'Left': __bLeft__, 'Right': __bRight__}

        #------- mesh examples ---------------------------------------------------------------------
        DICT = {
            'test0': boundaries_test0,
            'test2': boundaries_test0,
            'rand0': boundaries_test0,
        }
        #------- stK meshes ------------------------------------------------------------------------
        stK = [2, 4, 6, 8, 10, 12, 14, 16, 32, 48, 64 ,128]
        for key in stK:
            DICT['st' + str(key)] = boundaries_test0
        #------- stKd meshes -----------------------------------------------------------------------
        DICT.update(
            {
                'sqK8d2': boundaries_test1,
                'sqK16d2': boundaries_test1,
                'sqK32d2': boundaries_test1,
                'sqK48d2': boundaries_test1,
                'sqK64d2': boundaries_test1,
                'sqK128d2': boundaries_test1,
            }
        )

        return DICT




if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
