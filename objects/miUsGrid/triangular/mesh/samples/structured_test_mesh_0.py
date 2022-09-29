# -*- coding: utf-8 -*-
"""
With this script, we mesh a vtu file defining a structured triangular mesh for testing purposes.

@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/26 5:43 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')
import numpy as np
from root.config.main import rAnk, mAster_rank
from pyevtk.hl import unstructuredGridToVTK

def Generating_Structured_Triangular_Mesh(K):
    """The domain is [0,1]^2. It is divided into K*K uniform squares. And each square is then
    divided into two triangles.

    h = 1/K

    Parameters
    ----------
    K

    Returns
    -------

    """
    if rAnk == mAster_rank:
        x = np.tile(np.linspace(0, 1, K+1), (K+1))
        y = np.repeat(np.linspace(0, 1, K+1), (K+1))
        z = np.zeros((K+1)**2)
        points_matrix = np.arange(0, (K+1)**2).reshape([(K+1), (K+1)], order='F')

        CON = np.zeros((2*K**2, 3), dtype=int)

        for j in range(K):
            for i in range(K):
                Square = i + j * K

                points4 = (
                    points_matrix[i  , j  ],
                    points_matrix[i+1, j  ],
                    points_matrix[i  , j+1],
                    points_matrix[i+1, j+1]
                )

                c0 = Square * 2
                c1 = Square * 2 + 1

                CON[c0, :] = [points4[2], points4[0], points4[3]]
                CON[c1, :] = [points4[1], points4[3], points4[0]]

        CON = CON.ravel('C')
        offsets = np.arange(3, 3*2*K**2+3, 3)
        types = np.ones(2*K**2) * 5  # all triangles

        unstructuredGridToVTK(f'./objects/miUsGrid/triangular/mesh/samples/structured_tests/structured_test_mesh_01_K{K}',
                              x, y, z, CON, offsets, types)

if __name__ == "__main__":
    # mpiexec -n 4 python objects/miUsGrid/triangular/mesh/samples/structured_test_mesh_0.py
    Generating_Structured_Triangular_Mesh(2)