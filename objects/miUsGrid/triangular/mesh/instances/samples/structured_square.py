# -*- coding: utf-8 -*-
"""
With this script, we mesh a vtu file defining a structured triangular mesh for testing purposes.

@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/26 5:43 PM
"""
import os
import sys

if './' not in sys.path: sys.path.append('./')
import numpy as np
from root.config.main import RANK, MASTER_RANK
from pyevtk.hl import unstructuredGridToVTK



def Generating_Structured_Square_Mesh(K, d=2):
    """The domain is [-d/2,d/2]^2. It is divided into K*K uniform squares. And each square is then
    divided into two triangles.

    Parameters
    ----------
    K
    d

    Returns
    -------

    """
    if RANK == MASTER_RANK:

        x = np.tile(np.linspace(-d/2, d/2, K+1), (K+1))
        y = np.repeat(np.linspace(-d/2, d/2, K+1), (K+1))
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
        current_dir = os.path.dirname(__file__)

        unstructuredGridToVTK(current_dir + f'/structured_square/K{K}d{d}',
                              x, y, z, CON, offsets, types)



if __name__ == "__main__":
    # mpiexec -n 4 python objects/miUsGrid/triangular/mesh/instances/samples/structured_square.py
    Generating_Structured_Square_Mesh(48)