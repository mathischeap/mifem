# -*- coding: utf-8 -*-
"""
Yi Zhang
zhangyi_aero@hotmail.com
created at: 2/22/2023 12:41 PM
"""
import sys

if './' not in sys.path:
    sys.path.append('./')

import numpy as np

from __init__ import cscg3
from __init__ import components


def a(t, x, y, z): return 2 * t * np.sin(np.pi * x) * np.cos(1.111 * np.pi * y) * np.cos(np.pi * z)


def b(t, x, y, z): return t * np.sin(1.88 * np.pi * x) * np.cos(np.pi * y) * np.sin(0.998 * np.pi * z)


def c(t, x, y, z): return t * np.sin(2.09 * np.pi * x) * np.sin(0.567 * np.pi * y) * np.cos(np.pi * z)


def d(t, x, y, z): return 2 * t * np.cos(0.396 * np.pi * x) * np.cos(np.pi * y) * np.cos(np.pi * z)


def e(t, x, y, z): return t * np.cos(2.102 * np.pi * x) * np.cos(np.pi * y) * np.sin(np.pi * z)


def f(t, x, y, z): return t * np.cos(1.5 * np.pi * x) * np.cos(np.pi * y) * np.cos(1.35 * np.pi * z)


uu = components.ft3dw.vector(a, b, c)
BB = components.ft3dw.vector(d, e, f)

uXB = uu.cross_product(BB)


def cross_product_test():
    """"""

    mesh = cscg3.mesh('crazy', bounds=[(-1, 1) for _ in range(3)])([3, 3, 3], EDM='debug')
    space = cscg3.space('polynomials')([2,2,2])
    FC = cscg3.form(mesh, space)


    u = FC('2-f', hybrid=False)
    B = FC('2-f', hybrid=False)
    b = FC('2-f', hybrid=False)

    u.special.inner_curl__cross_product_2f(B, b)




if __name__ == '__main__':
    # mpiexec -n 4 python tests/objects/CSCG/_3d/unittests/standard_forms/inner_curl_X_2f_2f_1f.py

    cross_product_test()
