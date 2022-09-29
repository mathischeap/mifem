# -*- coding: utf-8 -*-
"""
Here we have constructed a random triangular mesh in the domain [0,1]^2 for testing purpose.

mpiexec -n 4 python objects/miUsGrid/triangular/__test__/Random/test_mesh.py

@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/20 3:45 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

import numpy as np

from root.config.main import rAnk, mAster_rank

from screws.miscellaneous.mios import remove
from pyevtk.hl import unstructuredGridToVTK
from objects.miUsGrid.triangular.mesh.main import miUsGrid_TriangularMesh

x = np.array([0, 0.3, 0.7, 1, 1, 0.5, 0, 0, 0.25, 0.32, 0.55, 1])
y = np.array([0, 0, 0, 0, 0.55, 0.44, 0.52, 1, 0.75, 1, 1, 1])
z = np.zeros(12)

CON = np.array(
    [0, 1, 6,
     5, 6, 1,
     1, 5, 2,
     2, 5, 4,
     3, 4, 2,
     5, 11, 4,
     11, 5, 10,
     10, 8, 5,
     8, 6, 5,
     8, 6, 7,
     7, 9, 8,
     9, 10, 8
     ]
)
offsets = np.array([3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36])
types = np.ones(12) * 5  # all triangles

if rAnk == mAster_rank:  # make sure only the master rank makes the vtu file.
    unstructuredGridToVTK('triangles', x, y, z, CON, offsets, types)

def bUpper(x, y): return x == 0
def bDown(x, y): return x == 1
def bLeft(x, y): return y == 0
def bRight(x, y): return y == 1

boundaries = {'Upper': bUpper, 'Down': bDown, 'Left': bLeft, 'Right':bRight}

mesh = miUsGrid_TriangularMesh('triangles.vtu', boundaries)

remove('triangles.vtu')