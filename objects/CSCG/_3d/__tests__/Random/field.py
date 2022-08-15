# -*- coding: utf-8 -*-
"""
INTRO

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
import sys
if './' not in sys.path: sys.path.append('./')
import random

from root.config.main import rAnk, mAster_rank, np, cOmm
from objects.CSCG._3d.fields.scalar.main import _3dCSCG_ScalarField
from objects.CSCG._3d.fields.vector.main import _3dCSCG_VectorField

from screws.miscellaneous.random_string.digits import randomStringDigits

def random_vector(mesh):
    """"""
    if rAnk == mAster_rank:
        a, b, c, d, e, f, g, h, i, j, k, l = [random.uniform(0.5, 1) for _ in range(12)]
        name = randomStringDigits(8)
    else:
        a, b, c, d, e, f, g, h, i, j, k, l = [None for _ in range(12)]
        name = None

    a, b, c, d, e, f, g, h, i, j, k, l = cOmm.bcast([a, b, c, d, e, f, g, h, i, j, k, l], root=mAster_rank)
    name = cOmm.bcast(name, root=mAster_rank)

    def u0(t, x, y, z): return np.sin(a*np.pi*x) * np.sin(b*np.pi*y) * np.cos(c*np.pi*z) + j*t
    def u1(t, x, y, z): return np.cos(d*np.pi*x) * np.sin(e*np.pi*y) * np.sin(f*np.pi*z) + k*t
    def u2(t, x, y, z): return np.sin(g*np.pi*x) * np.cos(h*np.pi*y) * np.sin(i*np.pi*z) + l*t

    return _3dCSCG_VectorField(mesh, (u0, u1, u2), name=name)

def random_scalar(mesh):
    """"""
    if rAnk == mAster_rank:
        a, b, c, d = [random.uniform(0.5, 1) for _ in range(4)]
        name = randomStringDigits(8)
    else:
        a, b, c, d = [None for _ in range(4)]
        name = None

    a, b, c, d = cOmm.bcast([a, b, c, d], root=mAster_rank)
    name = cOmm.bcast(name, root=mAster_rank)

    def p(t, x, y, z): return np.sin(a*np.pi*x) * np.sin(b*np.pi*y) * np.sin(d*np.pi*z)+ c*t

    return _3dCSCG_ScalarField(mesh, p, name=name)





if __name__ == '__main__':
    # mpiexec -n 4 python objects/CSCG/_3d/__tests__/Random/field.py
    random_vector(1)
