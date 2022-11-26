# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 11/21/2022 4:56 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from root.config.main import RANK, MASTER_RANK

from __init__ import components
ft3dw = components.ft3dw

from numpy import sin, cos, pi, exp
import numpy as np

def u(t, x, y, z):
    return sin(pi*x) * cos(pi*y) * cos(pi*z) * exp(2*t)
def ut(t, x, y, z):
    return sin(pi*x) * cos(pi*y) * cos(pi*z) * exp(2*t) * 2
def ux(t, x, y, z):
    return pi * cos(pi*x) * cos(pi*y) * cos(pi*z) * exp(2*t)
def uy(t, x, y, z):
    return - pi * sin(pi*x) * sin(pi*y) * cos(pi*z) * exp(2*t)
def uz(t, x, y, z):
    return - pi * sin(pi*x) * cos(pi*y) * sin(pi*z) * exp(2*t)

def v(t, x, y, z):
    return cos(pi*x) * sin(pi*y) * cos(pi*z) * exp(2*t)
def vt(t, x, y, z):
    return cos(pi*x) * sin(pi*y) * cos(pi*z) * exp(2*t) * 2
def vx(t, x, y, z):
    return - pi * sin(pi*x) * sin(pi*y) * cos(pi*z) * exp(2*t)
def vy(t, x, y, z):
    return pi * cos(pi*x) * cos(pi*y) * cos(pi*z) * exp(2*t)
def vz(t, x, y, z):
    return - pi * cos(pi*x) * sin(pi*y) * sin(pi*z) * exp(2*t)

def w(t, x, y, z):
    return cos(pi*x) * cos(pi*y) * sin(pi*z) * exp(2*t)
def wt(t, x, y, z):
    return cos(pi*x) * cos(pi*y) * sin(pi*z) * exp(2*t) * 2
def wx(t, x, y, z):
    return - pi * sin(pi*x) * cos(pi*y) * sin(pi*z) * exp(2*t)
def wy(t, x, y, z):
    return - pi * cos(pi*x) * sin(pi*y) * sin(pi*z) * exp(2*t)
def wz(t, x, y, z):
    return pi * cos(pi*x) * cos(pi*y) * cos(pi*z) * exp(2*t)

def f(t, x, y, z):
    return sin(pi*x) * sin(pi*y) * sin(pi*z) * exp(2*t)
def ft(t, x, y, z):
    return sin(pi*x) * sin(pi*y) * sin(pi*z) * exp(2*t) * 2
def fx(t, x, y, z):
    return pi * cos(pi*x) * sin(pi*y) * sin(pi*z) * exp(2*t)
def fy(t, x, y, z):
    return pi * sin(pi*x) * cos(pi*y) * sin(pi*z) * exp(2*t)
def fz(t, x, y, z):
    return pi * sin(pi*x) * sin(pi*y) * cos(pi*z) * exp(2*t)


def test_functions_time_plus_3d_wrappers_AKA_ft3dw():
    """"""
    if RANK == MASTER_RANK:
        print("ft2dw [test_functions_time_plus_3d_wrappers_AKA_ft3dw] ...... ", flush=True)

    t = np.random.rand(1)[0]
    x = np.random.rand(8)
    y = np.random.rand(8)
    z = np.random.rand(8)

    scalar = ft3dw.scalar(f)
    vector = ft3dw.vector(u, v, w)

    U, V, W, F = u(t,x,y,z), v(t,x,y,z), w(t,x,y,z), f(t,x,y,z)
    FX, FY, FZ = fx(t,x,y,z), fy(t,x,y,z), fz(t,x,y,z)
    UT, UX, UY, UZ = ut(t,x,y,z), ux(t,x,y,z), uy(t,x,y,z), uz(t,x,y,z)
    VT, VX, VY, VZ = vt(t,x,y,z), vx(t,x,y,z), vy(t,x,y,z), vz(t,x,y,z)
    WT, WX, WY, WZ = wt(t,x,y,z), wx(t,x,y,z), wy(t,x,y,z), wz(t,x,y,z)

    #------- test scalar -------------------------------------------------------------------
    SCALAR = scalar(t,x,y,z)
    np.testing.assert_array_almost_equal(SCALAR, F)

    sT = scalar.time_derivative
    np.testing.assert_array_almost_equal(sT(t,x,y,z), ft(t,x,y,z))
    sG = scalar.gradient
    sGv = sG(t,x,y,z)
    np.testing.assert_array_almost_equal(sGv[0], FX)
    np.testing.assert_array_almost_equal(sGv[1], FY)
    np.testing.assert_array_almost_equal(sGv[2], FZ)

    CONVECTION = scalar.convection_by(vector)(t,x,y,z)
    np.testing.assert_array_almost_equal(CONVECTION, U*FX + V*FY + W*FZ)

    sADD = sT + scalar
    sSub = sT - scalar
    sNeg = - sT
    ST = sT(t,x,y,z)
    np.testing.assert_array_almost_equal(sADD(t,x,y,z), ST + SCALAR)
    np.testing.assert_array_almost_equal(sSub(t,x,y,z), ST - SCALAR)
    np.testing.assert_array_almost_equal(sNeg(t,x,y,z), - ST)

    MULTIPLY = sT * scalar
    np.testing.assert_array_almost_equal(MULTIPLY(t,x,y,z), ST * SCALAR)

    #------------test vector ------------------------------------------------------------------
    VECTOR = vector(t,x,y,z)
    np.testing.assert_array_almost_equal(VECTOR[0], U)
    np.testing.assert_array_almost_equal(VECTOR[1], V)
    np.testing.assert_array_almost_equal(VECTOR[2], W)

    vDiv = vector.divergence
    vCurl = vector.curl
    vT = vector.time_derivative

    DIV = vDiv(t,x,y,z)
    CURL = vCurl(t,x,y,z)
    vTime = vT(t,x,y,z)
    np.testing.assert_array_almost_equal(DIV, UX + VY + WZ)
    np.testing.assert_array_almost_equal(CURL[0], WY-VZ)
    np.testing.assert_array_almost_equal(CURL[1], UZ-WX)
    np.testing.assert_array_almost_equal(CURL[2], VX-UY)
    np.testing.assert_array_almost_equal(vTime[0], UT)
    np.testing.assert_array_almost_equal(vTime[1], VT)
    np.testing.assert_array_almost_equal(vTime[2], WT)

    CONVECTION = vector.convection_by(vector)(t,x,y,z)
    np.testing.assert_array_almost_equal(CONVECTION[0], U*UX + V*UY + W*UZ)
    np.testing.assert_array_almost_equal(CONVECTION[1], U*VX + V*VY + W*VZ)
    np.testing.assert_array_almost_equal(CONVECTION[2], U*WX + V*WY + W*WZ)

    vADD = (vector + vT)(t,x,y,z)
    vSUB = (vector - vT)(t,x,y,z)
    vNEG = (- vT)(t,x,y,z)
    vDOT = vector.dot(vT)(t,x,y,z)
    np.testing.assert_array_almost_equal(vADD[0], U + vTime[0])
    np.testing.assert_array_almost_equal(vADD[1], V + vTime[1])
    np.testing.assert_array_almost_equal(vADD[2], W + vTime[2])
    np.testing.assert_array_almost_equal(vSUB[0], U - vTime[0])
    np.testing.assert_array_almost_equal(vSUB[1], V - vTime[1])
    np.testing.assert_array_almost_equal(vSUB[2], W - vTime[2])
    np.testing.assert_array_almost_equal(vNEG[0], - vTime[0])
    np.testing.assert_array_almost_equal(vNEG[1], - vTime[1])
    np.testing.assert_array_almost_equal(vNEG[2], - vTime[2])

    np.testing.assert_array_almost_equal(vDOT, U*vTime[0] + V*vTime[1] + W*vTime[2])

    return 1


if __name__ == '__main__':
    # mpiexec -n 4 python tests/components/ft3dw.py
    test_functions_time_plus_3d_wrappers_AKA_ft3dw()
