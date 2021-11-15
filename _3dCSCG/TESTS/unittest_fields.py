
"""
For testing fields

"""

import sys
if './' not in sys.path: sys.path.append('./')

from root.config import *
from SCREWS.miscellaneous import MyTimer
import random
from _3dCSCG.TESTS.random_objects import random_3D_mesh_and_space_of_total_load_around
from _3dCSCG.main import FormCaller, MeshGenerator, SpaceInvoker



def test_Form_NO0_3dCSCG_Field_numerical():
    """"""
    if rAnk == mAster_rank:
        load = random.randint(100, 499)
        print(f"-N- [test_Form_NO0_3dCSCG_Field_numerical] start at {MyTimer.current_time()} with load={load}.", flush=True)
    else:
        load= None
    load = cOmm.bcast(load, root=mAster_rank)
    mesh, space = random_3D_mesh_and_space_of_total_load_around(load, exclude_periodic=False)
    FC = FormCaller(mesh, space)

    t = random.random() * 10
    I, J, K = random.randint(2,10), random.randint(4,8), random.randint(3,9)
    x = np.linspace(-0.9-random.random()/10, 0.9+random.random()/10, I)
    y = np.linspace(-0.8-random.random()/5, 0.8+random.random()/5, J)
    z = np.linspace(-0.95-random.random()/20, 0.95+random.random()/20, K)



    #================ test for scalar ======================================================================
    def func(t, x, y, z): return np.sin(2*np.pi*x) * np.sin(np.pi*y) * np.sin(np.pi*z) * t
    SS = FC('scalar', func)
    def APt(t, x, y, z): return np.sin(2*np.pi*x) * np.sin(np.pi*y) * np.sin(np.pi*z) + 0*t
    def APx(t, x, y, z): return 2*np.pi*np.cos(2*np.pi*x) * np.sin(np.pi*y) * np.sin(np.pi*z) * t
    def APy(t, x, y, z): return np.pi*np.sin(2*np.pi*x) * np.cos(np.pi*y) * np.sin(np.pi*z) * t
    def APz(t, x, y, z): return np.pi*np.sin(2*np.pi*x) * np.sin(np.pi*y) * np.cos(np.pi*z) * t

    #---------- scalar gradient ---------
    numerical_gradient = SS.numerical.gradient
    numerical_gradient.current_time = t
    # in each core, current_time is different. This is bad for applications, but for tests, it is good.
    R_xyz, R_v = numerical_gradient.reconstruct(x, y, z)

    for i in R_xyz:
        xyz = R_xyz[i]
        exact_vx = APx(t, *xyz)
        exact_vy = APy(t, *xyz)
        exact_vz = APz(t, *xyz)
        numerical_vx, numerical_vy, numerical_vz = R_v[i]
        assert np.max(np.abs(exact_vx - numerical_vx)) < 1e-7, f"numerical partial x is not accurate enough."
        assert np.max(np.abs(exact_vy - numerical_vy)) < 1e-7, f"numerical partial y is not accurate enough."
        assert np.max(np.abs(exact_vz - numerical_vz)) < 1e-7, f"numerical partial z is not accurate enough."

    #---------- scalar time_derivative ---------
    time_derivative = SS.numerical.time_derivative
    time_derivative.current_time = t
    # in each core, current_time is different. This is bad for applications, but for tests, it is good.
    R_xyz, R_v = time_derivative.reconstruct(x, y, z)
    for i in R_xyz:
        xyz = R_xyz[i]
        exact_td = APt(t, *xyz)
        numerical_td = R_v[i]
        assert np.max(np.abs(exact_td - numerical_td)) < 1e-7, f"time derivative is not accurate enough."





    #================ test for vector ======================================================================
    def U(t, x, y, z): return 2 * t * np.sin(np.pi * x) * np.cos(np.pi * y) * np.cos(np.pi * z)
    def V(t, x, y, z): return t * np.sin(2 * np.pi * x) * np.cos(np.pi * y) * np.sin(np.pi * z)
    def W(t, x, y, z): return t * np.sin(2 * np.pi * x) * np.sin(np.pi * y) * np.cos(np.pi * z)

    SV = FC('vector', [U, V, W])

    def Ut(t, x, y, z): return 2 * np.sin(np.pi * x) * np.cos(np.pi * y) * np.cos(np.pi * z) + 0*t
    def Vt(t, x, y, z): return np.sin(2 * np.pi * x) * np.cos(np.pi * y) * np.sin(np.pi * z) + 0*t
    def Wt(t, x, y, z): return np.sin(2 * np.pi * x) * np.sin(np.pi * y) * np.cos(np.pi * z) + 0*t

    def Ux(t, x, y, z): return 2 * np.pi * t * np.cos(np.pi * x) * np.cos(np.pi * y) * np.cos(np.pi * z)
    def Uy(t, x, y, z): return - 2 * np.pi * t * np.sin(np.pi * x) * np.sin(np.pi * y) * np.cos(np.pi * z)
    def Uz(t, x, y, z): return - 2 * np.pi * t * np.sin(np.pi * x) * np.cos(np.pi * y) * np.sin(np.pi * z)

    def Vx(t, x, y, z): return 2 * np.pi * t * np.cos(2 * np.pi * x) * np.cos(np.pi * y) * np.sin(np.pi * z)
    def Vy(t, x, y, z): return - np.pi * t * np.sin(2 * np.pi * x) * np.sin(np.pi * y) * np.sin(np.pi * z)
    def Vz(t, x, y, z): return np.pi * t * np.sin(2 * np.pi * x) * np.cos(np.pi * y) * np.cos(np.pi * z)

    def Wx(t, x, y, z): return 2 * np.pi * t * np.cos(2 * np.pi * x) * np.sin(np.pi * y) * np.cos(np.pi * z)
    def Wy(t, x, y, z): return np.pi * t * np.sin(2 * np.pi * x) * np.cos(np.pi * y) * np.cos(np.pi * z)
    def Wz(t, x, y, z): return - np.pi * t * np.sin(2 * np.pi * x) * np.sin(np.pi * y) * np.sin(np.pi * z)

    #---------- vector time_derivative ---------
    time_derivative = SV.numerical.time_derivative
    time_derivative.current_time = t
    # in each core, current_time is different. This is bad for applications, but for tests, it is good.
    R_xyz, R_v = time_derivative.reconstruct(x, y, z)
    for i in R_xyz:
        xyz = R_xyz[i]
        exact_tdx = Ut(t, *xyz)
        exact_tdy = Vt(t, *xyz)
        exact_tdz = Wt(t, *xyz)
        numerical_tdx, numerical_tdy, numerical_tdz = R_v[i]
        assert np.max(np.abs(exact_tdx - numerical_tdx)) < 1e-7, f"time derivative x-component is not accurate enough."
        assert np.max(np.abs(exact_tdy - numerical_tdy)) < 1e-7, f"time derivative y-component is not accurate enough."
        assert np.max(np.abs(exact_tdz - numerical_tdz)) < 1e-7, f"time derivative z-component is not accurate enough."

    #---------- vector gradient ---------
    numerical_gradient = SV.numerical.gradient
    numerical_gradient.current_time = t
    # in each core, current_time is different. This is bad for applications, but for tests, it is good.
    R_xyz, R_v = numerical_gradient.reconstruct(x, y, z)

    for i in R_xyz:
        xyz = R_xyz[i]
        exact_Ux, exact_Uy, exact_Uz = Ux(t, *xyz), Uy(t, *xyz), Uz(t, *xyz)
        exact_Vx, exact_Vy, exact_Vz = Vx(t, *xyz), Vy(t, *xyz), Vz(t, *xyz)
        exact_Wx, exact_Wy, exact_Wz = Wx(t, *xyz), Wy(t, *xyz), Wz(t, *xyz)
        numerical_U, numerical_V, numerical_W = R_v[i]
        numerical_Ux, numerical_Uy, numerical_Uz = numerical_U
        numerical_Vx, numerical_Vy, numerical_Vz = numerical_V
        numerical_Wx, numerical_Wy, numerical_Wz = numerical_W

        assert np.max(np.abs(exact_Ux - numerical_Ux)) < 1e-7, f"numerical partial Ux is not accurate enough."
        assert np.max(np.abs(exact_Uy - numerical_Uy)) < 1e-7, f"numerical partial Uy is not accurate enough."
        assert np.max(np.abs(exact_Uz - numerical_Uz)) < 1e-7, f"numerical partial Uz is not accurate enough."
        assert np.max(np.abs(exact_Vx - numerical_Vx)) < 1e-7, f"numerical partial Vx is not accurate enough."
        assert np.max(np.abs(exact_Vy - numerical_Vy)) < 1e-7, f"numerical partial Vy is not accurate enough."
        assert np.max(np.abs(exact_Vz - numerical_Vz)) < 1e-7, f"numerical partial Vz is not accurate enough."
        assert np.max(np.abs(exact_Wx - numerical_Wx)) < 1e-7, f"numerical partial Wx is not accurate enough."
        assert np.max(np.abs(exact_Wy - numerical_Wy)) < 1e-7, f"numerical partial Wy is not accurate enough."
        assert np.max(np.abs(exact_Wz - numerical_Wz)) < 1e-7, f"numerical partial Wz is not accurate enough."

    #---------- vector curl ---------
    numerical_curl = SV.numerical.curl
    numerical_curl.current_time = t
    # in each core, current_time is different. This is bad for applications, but for tests, it is good.
    R_xyz, R_v = numerical_curl.reconstruct(x, y, z)

    for i in R_xyz:
        xyz = R_xyz[i]
        Cx = Wy(t, *xyz) - Vz(t, *xyz)
        Cy = Uz(t, *xyz) - Wx(t, *xyz)
        Cz = Vx(t, *xyz) - Uy(t, *xyz)
        numerical_Cx, numerical_Cy, numerical_Cz = R_v[i]
        assert np.max(np.abs(Cx - numerical_Cx)) < 1e-7, f"numerical x-component of curl is not accurate enough."
        assert np.max(np.abs(Cy - numerical_Cy)) < 1e-7, f"numerical y-component of curl is not accurate enough."
        assert np.max(np.abs(Cz - numerical_Cz)) < 1e-7, f"numerical z-component of curl is not accurate enough."

    #---------- vector divergence ---------
    numerical_divergence = SV.numerical.divergence
    numerical_divergence.current_time = t
    # in each core, current_time is different. This is bad for applications, but for tests, it is good.
    R_xyz, R_v = numerical_divergence.reconstruct(x, y, z)

    for i in R_xyz:
        xyz = R_xyz[i]
        D = Ux(t, *xyz) + Vy(t, *xyz) + Wz(t, *xyz)
        numerical_D = R_v[i][0]
        assert np.max(np.abs(D - numerical_D)) < 1e-7, f"numerical divergence is not accurate enough."





    #================ test for tensor ======================================================================
    def T00(t, x, y, z): return 2 * t * np.sin(np.pi * x) * np.cos(np.pi * y) * np.cos(np.pi * z)
    def T01(t, x, y, z): return t * np.sin(2 * np.pi * x) * np.cos(np.pi * y) * np.sin(np.pi * z)
    def T02(t, x, y, z): return t * np.sin(2 * np.pi * x) * np.sin(np.pi * y) * np.cos(np.pi * z)

    def T10(t, x, y, z): return 2 * t * np.cos(np.pi * x) * np.cos(np.pi * y) * np.cos(np.pi * z)
    def T11(t, x, y, z): return t * np.cos(2 * np.pi * x) * np.cos(np.pi * y) * np.sin(np.pi * z)
    def T12(t, x, y, z): return t * np.cos(2 * np.pi * x) * np.sin(np.pi * y) * np.cos(np.pi * z)

    def T20(t, x, y, z): return t * np.sin(np.pi * x) * np.cos(np.pi * y) * np.sin(np.pi * z)
    def T21(t, x, y, z): return t * np.sin(2 * np.pi * x) * np.cos(np.pi * y) * np.sin(np.pi * z)
    def T22(t, x, y, z): return t * np.sin(2 * np.pi * x) * np.sin(np.pi * y) * np.sin(np.pi * z)

    ST = FC('tensor', ([T00, T01, T02], [T10, T11, T12], [T20, T21, T22]))

    #---------- tensor time_derivative ---------
    def T00t(t, x, y, z): return 2 * np.sin(np.pi * x) * np.cos(np.pi * y) * np.cos(np.pi * z) + 0 * t
    def T01t(t, x, y, z): return np.sin(2 * np.pi * x) * np.cos(np.pi * y) * np.sin(np.pi * z) + 0 * t
    def T02t(t, x, y, z): return np.sin(2 * np.pi * x) * np.sin(np.pi * y) * np.cos(np.pi * z) + 0 * t

    def T10t(t, x, y, z): return 2 * np.cos(np.pi * x) * np.cos(np.pi * y) * np.cos(np.pi * z) + 0 * t
    def T11t(t, x, y, z): return np.cos(2 * np.pi * x) * np.cos(np.pi * y) * np.sin(np.pi * z) + 0 * t
    def T12t(t, x, y, z): return np.cos(2 * np.pi * x) * np.sin(np.pi * y) * np.cos(np.pi * z) + 0 * t

    def T20t(t, x, y, z): return np.sin(np.pi * x) * np.cos(np.pi * y) * np.sin(np.pi * z) + 0 * t
    def T21t(t, x, y, z): return np.sin(2 * np.pi * x) * np.cos(np.pi * y) * np.sin(np.pi * z) + 0 * t
    def T22t(t, x, y, z): return np.sin(2 * np.pi * x) * np.sin(np.pi * y) * np.sin(np.pi * z) + 0 * t
    time_derivative = ST.numerical.time_derivative
    time_derivative.current_time = t
    # in each core, current_time is different. This is bad for applications, but for tests, it is good.
    R_xyz, R_v = time_derivative.reconstruct(x, y, z)
    for i in R_xyz:
        xyz = R_xyz[i]
        AT00, AT01, AT02 = T00t(t, *xyz), T01t(t, *xyz), T02t(t, *xyz)
        AT10, AT11, AT12 = T10t(t, *xyz), T11t(t, *xyz), T12t(t, *xyz)
        AT20, AT21, AT22 = T20t(t, *xyz), T21t(t, *xyz), T22t(t, *xyz)
        numerical_T0, numerical_T1, numerical_T2 = R_v[i]
        numerical_T00, numerical_T01, numerical_T02 = numerical_T0
        numerical_T10, numerical_T11, numerical_T12 = numerical_T1
        numerical_T20, numerical_T21, numerical_T22 = numerical_T2
        assert np.max(np.abs(AT00 - numerical_T00)) < 1e-7, f"time derivative 00-component is not accurate enough."
        assert np.max(np.abs(AT01 - numerical_T01)) < 1e-7, f"time derivative 01-component is not accurate enough."
        assert np.max(np.abs(AT02 - numerical_T02)) < 1e-7, f"time derivative 02-component is not accurate enough."
        assert np.max(np.abs(AT10 - numerical_T10)) < 1e-7, f"time derivative 10-component is not accurate enough."
        assert np.max(np.abs(AT11 - numerical_T11)) < 1e-7, f"time derivative 11-component is not accurate enough."
        assert np.max(np.abs(AT12 - numerical_T12)) < 1e-7, f"time derivative 12-component is not accurate enough."
        assert np.max(np.abs(AT20 - numerical_T20)) < 1e-7, f"time derivative 20-component is not accurate enough."
        assert np.max(np.abs(AT21 - numerical_T21)) < 1e-7, f"time derivative 21-component is not accurate enough."
        assert np.max(np.abs(AT22 - numerical_T22)) < 1e-7, f"time derivative 22-component is not accurate enough."

    #---------- tensor divergence ---------
    def T00x(t, x, y, z): return np.pi * 2 * t * np.cos(np.pi * x) * np.cos(np.pi * y) * np.cos(np.pi * z)
    def T01y(t, x, y, z): return - np.pi * t * np.sin(2 * np.pi * x) * np.sin(np.pi * y) * np.sin(np.pi * z)
    def T02z(t, x, y, z): return - np.pi * t * np.sin(2 * np.pi * x) * np.sin(np.pi * y) * np.sin(np.pi * z)

    def T10x(t, x, y, z): return - np.pi * 2 * t * np.sin(np.pi * x) * np.cos(np.pi * y) * np.cos(np.pi * z)
    def T11y(t, x, y, z): return - np.pi * t * np.cos(2 * np.pi * x) * np.sin(np.pi * y) * np.sin(np.pi * z)
    def T12z(t, x, y, z): return - np.pi * t * np.cos(2 * np.pi * x) * np.sin(np.pi * y) * np.sin(np.pi * z)

    def T20x(t, x, y, z): return np.pi * t * np.cos(np.pi * x) * np.cos(np.pi * y) * np.sin(np.pi * z)
    def T21y(t, x, y, z): return - np.pi * t * np.sin(2 * np.pi * x) * np.sin(np.pi * y) * np.sin(np.pi * z)
    def T22z(t, x, y, z): return np.pi * t * np.sin(2 * np.pi * x) * np.sin(np.pi * y) * np.cos(np.pi * z)

    numerical_divergence = ST.numerical.divergence
    numerical_divergence.current_time = t
    # in each core, current_time is different. This is bad for applications, but for tests, it is good.
    R_xyz, R_v = numerical_divergence.reconstruct(x, y, z)

    for i in R_xyz:
        xyz = R_xyz[i]
        D0 = T00x(t, *xyz) + T01y(t, *xyz) + T02z(t, *xyz)
        D1 = T10x(t, *xyz) + T11y(t, *xyz) + T12z(t, *xyz)
        D2 = T20x(t, *xyz) + T21y(t, *xyz) + T22z(t, *xyz)
        numerical_D0, numerical_D1, numerical_D2 = R_v[i]
        assert np.max(np.abs(D0 - numerical_D0)) < 1e-7, f"numerical x-component divergence is not accurate enough."
        assert np.max(np.abs(D1 - numerical_D1)) < 1e-7, f"numerical y-component divergence is not accurate enough."
        assert np.max(np.abs(D2 - numerical_D2)) < 1e-7, f"numerical z-component divergence is not accurate enough."

    return 1


def test_Form_NO1_3dCSCG_VectorField():
    """"""
    if rAnk == mAster_rank:
        print(f"-V- [test_Form_NO1_3dCSCG_VectorField]...", flush=True)
    def w0(t, x, y, z): return -np.pi * np.sin(x) * np.cos(y) * np.cos(z) + t
    def w1(t, x, y, z): return -np.pi * np.cos(x) * np.sin(y) * np.cos(z) + t
    def w2(t, x, y, z): return -np.pi * np.cos(x) * np.sin(y) * np.sin(z) + t
    def u0(t, x, y, z): return np.sin(np.pi*x) * np.sin(y) * np.cos(np.pi*z) + t
    def u1(t, x, y, z): return np.cos(x) * np.sin(np.pi*y) * np.sin(np.pi*z) + t
    def u2(t, x, y, z): return np.sin(np.pi*x) * np.cos(np.pi*y) * np.sin(z) + t
    # x = w X u
    def x0(t, x, y, z): return w1(t, x, y, z) * u2(t, x, y, z) - w2(t, x, y, z) * u1(t, x, y, z)
    def x1(t, x, y, z): return w2(t, x, y, z) * u0(t, x, y, z) - w0(t, x, y, z) * u2(t, x, y, z)
    def x2(t, x, y, z): return w0(t, x, y, z) * u1(t, x, y, z) - w1(t, x, y, z) * u0(t, x, y, z)

    t = random.random() * 10
    I, J, K = random.randint(2,10), random.randint(4,8), random.randint(3,9)
    x = np.linspace(-0.9-random.random()/10, 0.9+random.random()/10, I)
    y = np.linspace(-0.8-random.random()/5, 0.8+random.random()/5, J)
    z = np.linspace(-0.9-random.random()/10, 0.9+random.random()/10, K)

    mesh = MeshGenerator('crazy', c=0.0, bounds=([0,1],[0,1],[0,1]))([2,3,1], EDM='debug')
    space = SpaceInvoker('polynomials')([('Lobatto', 2), ('Lobatto', 2), ('Lobatto', 3)])
    FC = FormCaller(mesh, space)

    W = FC('vector', (w0, w1, w2))
    U = FC('vector', (u0, u1, u2))

    # ---------- neg ----------------------------------------------------

    mW = - W
    mW.current_time = t
    R_xyz, R_v = mW.reconstruct(x, y, z)
    for i in R_xyz:
        xyz = R_xyz[i]
        Ax0, Ax1, Ax2 = -w0(t, *xyz), -w1(t, *xyz), -w2(t, *xyz)
        Cx0, Cx1, Cx2 = R_v[i]
        assert np.max(np.abs(Ax0 - Cx0)) < 1e-10, f"neg x-component is not accurate enough."
        assert np.max(np.abs(Ax1 - Cx1)) < 1e-10, f"neg y-component is not accurate enough."
        assert np.max(np.abs(Ax2 - Cx2)) < 1e-10, f"neg z-component is not accurate enough."

    # ----------- sub ------------------------------------------------
    S = W - U
    S.current_time = t
    R_xyz, R_v = S.reconstruct(x, y, z)
    for i in R_xyz:
        xyz = R_xyz[i]
        Ax0, Ax1, Ax2 = w0(t, *xyz) - u0(t, *xyz), w1(t, *xyz) - u1(t, *xyz), w2(t, *xyz) - u2(t, *xyz)
        Cx0, Cx1, Cx2 = R_v[i]
        assert np.max(np.abs(Ax0 - Cx0)) < 1e-10, f"sub x-component is not accurate enough."
        assert np.max(np.abs(Ax1 - Cx1)) < 1e-10, f"sub y-component is not accurate enough."
        assert np.max(np.abs(Ax2 - Cx2)) < 1e-10, f"sub z-component is not accurate enough."


    # ----------- add ------------------------------------------------
    A = W + U
    A.current_time = t
    R_xyz, R_v = A.reconstruct(x, y, z)
    for i in R_xyz:
        xyz = R_xyz[i]
        Ax0, Ax1, Ax2 = w0(t, *xyz) + u0(t, *xyz), w1(t, *xyz) + u1(t, *xyz), w2(t, *xyz) + u2(t, *xyz)
        Cx0, Cx1, Cx2 = R_v[i]
        assert np.max(np.abs(Ax0 - Cx0)) < 1e-10, f"add x-component is not accurate enough."
        assert np.max(np.abs(Ax1 - Cx1)) < 1e-10, f"add y-component is not accurate enough."
        assert np.max(np.abs(Ax2 - Cx2)) < 1e-10, f"add z-component is not accurate enough."

    # ------ cross product ----------------------------------------
    X = W.DO.cross_product(U)
    X.current_time = t
    R_xyz, R_v = X.reconstruct(x, y, z)

    for i in R_xyz:
        xyz = R_xyz[i]
        Ax0, Ax1, Ax2 = x0(t, *xyz), x1(t, *xyz), x2(t, *xyz)
        Cx0, Cx1, Cx2 = R_v[i]
        assert np.max(np.abs(Ax0 - Cx0)) < 1e-10, f"cross product x-component is not accurate enough."
        assert np.max(np.abs(Ax1 - Cx1)) < 1e-10, f"cross product y-component is not accurate enough."
        assert np.max(np.abs(Ax2 - Cx2)) < 1e-10, f"cross product z-component is not accurate enough."



    return 1


def test_Form_NO2_3dCSCG_ScalarField():
    """"""
    if rAnk == mAster_rank:
        print(f"-S- [test_Form_NO2_3dCSCG_ScalarField]...", flush=True)

    def www(t, x, y, z): return -np.pi * np.sin(x) * np.cos(y) * np.cos(z) + np.cos(t)
    def uuu(t, x, y, z): return np.sin(2*np.pi*x) * np.sin(3*y) * np.cos(np.pi*0.125*z) + np.sin(t)


    t = random.random() * 10
    I, J, K = random.randint(2,10), random.randint(4,8), random.randint(3,9)
    x = np.linspace(-0.9-random.random()/10, 0.9+random.random()/10, I)
    y = np.linspace(-0.8-random.random()/5, 0.8+random.random()/5, J)
    z = np.linspace(-0.9-random.random()/10, 0.9+random.random()/10, K)

    mesh = MeshGenerator('crazy', c=0.0, bounds=([0,1],[0,1],[0,1]))([2,3,1], EDM='debug')
    space = SpaceInvoker('polynomials')([('Lobatto', 2), ('Lobatto', 2), ('Lobatto', 3)])
    FC = FormCaller(mesh, space)

    A = FC('scalar', www)
    B = FC('scalar', uuu)

    #---- neg --------------------------------------------------------------------
    X = -A
    X.current_time = t
    R_xyz, R_v = X.reconstruct(x, y, z)
    for i in R_xyz:
        xyz = R_xyz[i]
        Ax = -www(t, *xyz)
        Cx = R_v[i][0]
        assert np.max(np.abs(Ax - Cx)) < 1e-10, f"neg is not accurate enough."

    #---- sub --------------------------------------------------------------------
    X = A - B
    X.current_time = t
    R_xyz, R_v = X.reconstruct(x, y, z)
    for i in R_xyz:
        xyz = R_xyz[i]
        Ax = www(t, *xyz) - uuu(t, *xyz)
        Cx = R_v[i][0]
        assert np.max(np.abs(Ax - Cx)) < 1e-10, f"sub is not accurate enough."


    #---- add --------------------------------------------------------------------
    X = A + B
    X.current_time = t
    R_xyz, R_v = X.reconstruct(x, y, z)
    for i in R_xyz:
        xyz = R_xyz[i]
        Ax = www(t, *xyz) + uuu(t, *xyz)
        Cx = R_v[i][0]
        assert np.max(np.abs(Ax - Cx)) < 1e-10, f"add is not accurate enough."

    return 1



if __name__ == '__main__':
    # mpiexec -n 6 python _3dCSCG\TESTS\unittest_fields.py
    test_Form_NO1_3dCSCG_VectorField()
    test_Form_NO0_3dCSCG_Field_numerical()
    test_Form_NO2_3dCSCG_ScalarField()