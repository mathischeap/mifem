# -*- coding: utf-8 -*-
"""
mpiexec -n 4 python objects/miUsGrid/triangular/__test__/Random/test_mesh_2.py

@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 10/2/2022 4:24 PM
"""
import sys

if './' not in sys.path:
    sys.path.append('/')
import numpy as np

from root.config.main import RANK, MASTER_RANK

from components.miscellaneous.mios import remove
from pyevtk.hl import unstructuredGridToVTK
from objects.miUsGrid.triangular.mesh.main import miUsGrid_TriangularMesh

x = np.array([0, 1, 1, 0])
y = np.array([0, 0, 1, 1])
z = np.zeros(4, dtype=int)

CON = np.array(
    [0, 1, 2,
     0, 2, 3
     ]
)
offsets = np.array([3, 6])
types = np.ones(2) * 5  # all triangles
if RANK == MASTER_RANK:  # make sure only the master rank makes the vtu file.
    unstructuredGridToVTK('test2', x, y, z, CON, offsets, types)


# noinspection PyUnusedLocal
def bUpper(x, y): return x == 0


# noinspection PyUnusedLocal
def bDown(x, y): return x == 1


# noinspection PyUnusedLocal
def bLeft(x, y): return y == 0


# noinspection PyUnusedLocal
def bRight(x, y): return y == 1


boundaries = {'Upper': bUpper, 'Down': bDown, 'Left': bLeft, 'Right': bRight}

mesh = miUsGrid_TriangularMesh('test2.vtu', boundaries)

remove('test2.vtu')
