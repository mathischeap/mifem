# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/8/10 23:53
"""
import sys
if './' not in sys.path:
    sys.path.append('/')
import random

from root.config.main import RANK, MASTER_RANK, np, COMM
from objects.CSCG._2d.fields.scalar.main import _2dCSCG_ScalarField
from objects.CSCG._2d.fields.vector.main import _2dCSCG_VectorField

from components.miscellaneous.randomString.digits import randomStringDigits


def random_vector(mesh):
    """"""
    if RANK == MASTER_RANK:
        a, b, j, k, d, e = [random.uniform(0.5, 1) for _ in range(6)]
        name = randomStringDigits(8)
    else:
        a, b, j, k, d, e = [None for _ in range(6)]
        name = None

    a, b, j, k, d, e = COMM.bcast([a, b, j, k, d, e], root=MASTER_RANK)
    name = COMM.bcast(name, root=MASTER_RANK)

    def u0(t, x, y): return np.sin(a*np.pi*x) * np.sin(b*np.pi*y) + j*t
    def u1(t, x, y): return np.cos(d*np.pi*x) * np.sin(e*np.pi*y) + k*t

    return _2dCSCG_VectorField(mesh, (u0, u1), name=name)


def random_scalar(mesh):
    """"""
    if RANK == MASTER_RANK:
        a, b, c = [random.uniform(0.5, 1) for _ in range(3)]
        name = randomStringDigits(8)
    else:
        a, b, c = [None for _ in range(3)]
        name = None

    a, b, c = COMM.bcast([a, b, c], root=MASTER_RANK)
    name = COMM.bcast(name, root=MASTER_RANK)

    def p(t, x, y): return np.sin(a*np.pi*x) * np.sin(b*np.pi*y) + c*t

    return _2dCSCG_ScalarField(mesh, p, name=name)


if __name__ == '__main__':
    # mpiexec -n 4 python objects/CSCG/_2d/__tests__/Random/field.py
    pass
