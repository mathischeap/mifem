# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 11/21/2022 3:55 PM
"""
import sys

if './' not in sys.path:
    sys.path.append('./')

from root.config.main import RANK, MASTER_RANK

from __init__ import components
ft2dw = components.ft2dw

from numpy import sin, cos, pi, exp
import numpy as np


def u(t, x, y):
    return cos(pi * x) * sin(pi * y) * exp(2*t)


def v(t, x, y):
    return sin(pi * x) * cos(pi * y) * exp(2*t)


def u_t(t, x, y):
    return cos(pi * x) * sin(pi * y) * exp(2*t) * 2


def v_t(t, x, y):
    return sin(pi * x) * cos(pi * y) * exp(2*t) * 2


def u_x(t, x, y):
    return - pi * sin(pi * x) * sin(pi * y) * exp(2*t)


def u_y(t, x, y):
    return pi * cos(pi * x) * cos(pi * y) * exp(2*t)


def v_x(t, x, y):
    return pi * cos(pi * x) * cos(pi * y) * exp(2*t)


def v_y(t, x, y):
    return - pi * sin(pi * x) * sin(pi * y) * exp(2*t)



def fs(t, x, y):
    return sin(pi * x) * sin(pi * y) * exp(2*t)


def fs_x(t, x, y):
    return pi * cos(pi * x) * sin(pi * y) * exp(2*t)


def fs_y(t, x, y):
    return pi * sin(pi * x) * cos(pi * y) * exp(2*t)


def fs_t(t, x, y):
    return sin(pi * x) * sin(pi * y) * 2 * exp(2*t)


def test_functions_time_plus_2d_wrappers_AKA_ft2dw():
    """"""
    if RANK == MASTER_RANK:
        print("ft2dw [test_functions_time_plus_2d_wrappers_AKA_ft2dw] ...... ", flush=True)

    x = np.random.rand(10)
    y = np.random.rand(10)
    t = np.random.rand(1)[0]

    scalar = ft2dw.scalar(fs)
    vector = ft2dw.vector(u, v)

    U_txy = u(t, x, y)
    V_txy = v(t, x, y)

    # -------- test scalar -------------------------------------------------------------------
    grad_scalar = scalar.gradient
    pt_scalar = scalar.time_derivative
    curl_scalar = scalar.curl

    SCALAR = scalar(t, x, y)
    np.testing.assert_array_almost_equal(SCALAR, fs(t, x, y))

    gS = grad_scalar(t, x, y)
    np.testing.assert_array_almost_equal(fs_x(t, x, y), gS[0])
    np.testing.assert_array_almost_equal(fs_y(t, x, y), gS[1])

    tS = pt_scalar(t, x, y)
    np.testing.assert_array_almost_equal(tS, fs_t(t, x, y))

    cS = curl_scalar(t, x, y)
    np.testing.assert_array_almost_equal(fs_y(t, x, y), cS[0])
    np.testing.assert_array_almost_equal(-fs_x(t, x, y), cS[1])

    sPlus = scalar + pt_scalar
    sSub = scalar - pt_scalar

    txy_plus = fs(t, x, y) + fs_t(t, x, y)
    np.testing.assert_array_almost_equal(sPlus(t, x, y), txy_plus)
    np.testing.assert_array_almost_equal(sSub(t, x, y), fs(t, x, y) - fs_t(t, x, y))

    negS = - scalar
    np.testing.assert_array_almost_equal(negS(t, x, y), - fs(t, x, y))

    MULTIPLY = sPlus * scalar
    np.testing.assert_array_almost_equal(MULTIPLY(t, x, y), fs(t, x, y) * txy_plus)

    convection = scalar.convection_by(vector)
    CON = convection(t, x, y)
    np.testing.assert_array_almost_equal(CON, U_txy * fs_x(t, x, y) + V_txy * fs_y(t, x, y))

    # -------- test vector -------------------------------------------------------------------

    V = vector(t, x, y)
    np.testing.assert_array_almost_equal(V[0], U_txy)
    np.testing.assert_array_almost_equal(V[1], V_txy)
    vT = vector.time_derivative
    VT = vT(t, x, y)
    np.testing.assert_array_almost_equal(VT[0], u_t(t, x, y))
    np.testing.assert_array_almost_equal(VT[1], v_t(t, x, y))

    div_vector = vector.divergence
    np.testing.assert_array_almost_equal(div_vector(t, x, y), u_x(t, x, y)+v_y(t, x, y))

    rot_vector = vector.rot
    np.testing.assert_array_almost_equal(rot_vector(t, x, y), v_x(t, x, y)-u_y(t, x, y))

    convection = vector.convection_by(vector)(t, x, y)
    np.testing.assert_array_almost_equal(convection[0], U_txy*u_x(t, x, y) + V_txy*u_y(t, x, y))
    np.testing.assert_array_almost_equal(convection[1], U_txy*v_x(t, x, y) + V_txy*v_y(t, x, y))

    vPlus = vector + vT
    v_Sub = vector - vT
    v_neg = - vT
    vP = vPlus(t, x, y)
    vS = v_Sub(t, x, y)
    vN = v_neg(t, x, y)
    np.testing.assert_array_almost_equal(vP[0], U_txy + VT[0])
    np.testing.assert_array_almost_equal(vP[1], V_txy + VT[1])
    np.testing.assert_array_almost_equal(vS[0], U_txy - VT[0])
    np.testing.assert_array_almost_equal(vS[1], V_txy - VT[1])
    np.testing.assert_array_almost_equal(vN[0], - VT[0])
    np.testing.assert_array_almost_equal(vN[1], - VT[1])

    DOT = vector.dot(vPlus)(t, x, y)
    np.testing.assert_array_almost_equal(DOT, (U_txy * vP[0]) + (V_txy * vP[1]))

    return 1


if __name__ == '__main__':
    # mpiexec -n 4 python tests/components/ft2dw.py
    test_functions_time_plus_2d_wrappers_AKA_ft2dw()
