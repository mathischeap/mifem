# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 10/2/2022 1:00 AM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from objects.miUsGrid.triangular.mesh.samples.rand_mesh_0 import generate_rand_mesh_0
from objects.miUsGrid.triangular.mesh.main import miUsGrid_TriangularMesh

generate_rand_mesh_0()

def bUpper(x, y): return x == 0
def bDown(x, y): return x == 1
def bLeft(x, y): return y == 0
def bRight(x, y): return y == 1

boundaries = {'Upper': bUpper, 'Down': bDown, 'Left': bLeft, 'Right':bRight}

mesh = miUsGrid_TriangularMesh(str(__file__).split('object')[0] +
        "objects/miUsGrid/triangular/mesh/samples/rand0.vtu", boundaries)
