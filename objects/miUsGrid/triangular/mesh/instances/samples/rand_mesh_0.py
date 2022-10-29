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

from root.config.main import RANK, MASTER_RANK, COMM
from random import uniform

from pyevtk.hl import unstructuredGridToVTK

def generate_rand_mesh_0():

    COMM.barrier()

    if RANK == MASTER_RANK:  # make sure only the master rank makes the vtu file.
        x = np.array([0, uniform(0.25, 0.35),
                      uniform(0.65, 0.75), 1, 1, uniform(0.4, 0.6), 0, 0,
                      uniform(0.22, 0.28), uniform(0.3, 0.4), uniform(0.5, 0.65), 1])
        y = np.array([0, 0, 0, 0, uniform(0.5, 0.6), uniform(0.4, 0.5),
                      uniform(0.5, 0.55), 1, uniform(0.7, 0.8), 1, 1, 1])
        z = np.zeros(12)

        CON = np.concatenate(
             np.random.permutation([
                 np.random.permutation([0, 1, 6]),
                 np.random.permutation([5, 6, 1]),
                 np.random.permutation([1, 5, 2]),
                 np.random.permutation([2, 5, 4]),
                 np.random.permutation([3, 4, 2]),
                 np.random.permutation([5, 11, 4]),
                 np.random.permutation([11, 5, 10]),
                 np.random.permutation([10, 8, 5]),
                 np.random.permutation([8, 6, 5]),
                 np.random.permutation([8, 6, 7]),
                 np.random.permutation([7, 9, 8]),
                 np.random.permutation([9, 10, 8])
                ])
        )
        offsets = np.array([3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36])
        types = np.ones(12) * 5  # all triangles
        base_path = str(__file__).split('rand_mesh_0')[0]
        unstructuredGridToVTK(base_path + 'rand0', x, y, z, CON, offsets, types)