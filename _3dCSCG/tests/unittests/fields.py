
"""
For testing fields

"""

import sys
if './' not in sys.path: sys.path.append('./')

from root.config.main import *
from screws.miscellaneous.timer import MyTimer
import random
from _3dCSCG.tests.random_objects.form_caller import random_3D_FormCaller_of_total_load_around



def test_Form_NO0_3dCSCG_Field_numerical():
    """"""
    if rAnk == mAster_rank:
        load = random.randint(100, 499)
        print(f"-N- [test_Form_NO0_3dCSCG_Field_numerical] {MyTimer.current_time()} with load={load}.", flush=True)
    else:
        load= None
    load = cOmm.bcast(load, root=mAster_rank)
    FC = random_3D_FormCaller_of_total_load_around(load, exclude_periodic=False)

    t = random.random() * 10
    I, J, K = random.randint(2,10), random.randint(4,8), random.randint(3,9)
    xi = np.linspace(-0.9-random.random()/10, 0.9+random.random()/10, I)
    et = np.linspace(-0.8-random.random()/5, 0.8+random.random()/5, J)
    sg = np.linspace(-0.95-random.random()/20, 0.95+random.random()/20, K)
    x, y, z = np.meshgrid(xi, et, sg, indexing='ij')

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
    def w0(t, x, y, z): return -np.pi * np.sin(x) * np.cos(y) * np.cos(z) + np.sin(t)
    def w1(t, x, y, z): return -np.pi * np.cos(x) * np.sin(y) * np.cos(z) * np.sin(2*t)
    def w2(t, x, y, z): return -np.pi * np.cos(x) * np.sin(y) * np.sin(z) / (1.5-np.cos(t/2))
    def u0(t, x, y, z): return np.sin(np.pi*x) * np.sin(y) * np.cos(np.pi*z) * 5 * np.sin(t)
    def u1(t, x, y, z): return np.cos(x) * np.sin(np.pi*y) * np.sin(np.pi*z) * 5 * np.cos(t)
    def u2(t, x, y, z): return np.sin(np.pi*x) * np.cos(np.pi*y) * np.sin(z) + 2 * np.cos(t)
    # x = w X u
    def x0(t, x, y, z): return w1(t, x, y, z) * u2(t, x, y, z) - w2(t, x, y, z) * u1(t, x, y, z)
    def x1(t, x, y, z): return w2(t, x, y, z) * u0(t, x, y, z) - w0(t, x, y, z) * u2(t, x, y, z)
    def x2(t, x, y, z): return w0(t, x, y, z) * u1(t, x, y, z) - w1(t, x, y, z) * u0(t, x, y, z)

    t = random.random() * 10
    t1 = random.random() * 10
    t2 = random.random() * 10
    t3 = random.random() * 10
    I, J, K = random.randint(2,10), random.randint(4,8), random.randint(3,9)
    xi = np.linspace(-0.9-random.random()/10, 0.9+random.random()/10, I)
    et = np.linspace(-0.8-random.random()/5, 0.8+random.random()/5, J)
    sg = np.linspace(-0.9-random.random()/10, 0.9+random.random()/10, K)
    x, y, z = np.meshgrid(xi, et, sg, indexing='ij')

    if rAnk == mAster_rank:
        load = random.randint(100, 499)
    else:
        load= None
    load = cOmm.bcast(load, root=mAster_rank)

    #----------------- use crazy mesh --------------------------------------------------------
    FC = random_3D_FormCaller_of_total_load_around(load, mesh_pool=('crazy',))
    W = FC('vector', (w0, w1, w2))
    W.current_time = t
    # ------ norm component ------------------------------------------------------------------
    N = W.components.norm
    N.current_time = t + t1
    assert N.current_time == W.current_time + t1, f"N and W has decoupled!"
    assert N.standard_properties.name == "norm-component-of-" + W.standard_properties.name, f"The naming rule."
    R_xyz, R_v = N.reconstruct(xi, et, sg, i='on_mesh_boundaries')
    for i in R_xyz:
        xyz = R_xyz[i]
        te = N.mesh.trace.elements[i]
        side = te.CHARACTERISTIC_side
        if side in 'NS':
            np.testing.assert_array_almost_equal(R_v[i][0] - w0(t+t1, *xyz), 0)
            np.testing.assert_array_almost_equal(R_v[i][1], 0)
            np.testing.assert_array_almost_equal(R_v[i][2], 0)
        elif side in 'WE':
            np.testing.assert_array_almost_equal(R_v[i][0], 0)
            np.testing.assert_array_almost_equal(R_v[i][1] - w1(t+t1, *xyz), 0)
            np.testing.assert_array_almost_equal(R_v[i][2], 0)
        elif side in 'BF':
            np.testing.assert_array_almost_equal(R_v[i][0], 0)
            np.testing.assert_array_almost_equal(R_v[i][1], 0)
            np.testing.assert_array_almost_equal(R_v[i][2] - w2(t+t1, *xyz), 0)
        else:
            raise Exception()

    # T_para component ---------------------------------------------------------------------------
    T_para = W.components.T_para
    T_para.current_time = t + t2
    assert T_para.current_time == W.current_time + t2, f"T_para and W has decoupled!"
    assert T_para.standard_properties.name == "T-para-component-of-" + W.standard_properties.name, f"The naming rule."
    R_xyz, R_v = T_para.reconstruct(xi, et, sg, i='on_mesh_boundaries')
    for i in R_xyz:
        xyz = R_xyz[i]
        te = T_para.mesh.trace.elements[i]
        side = te.CHARACTERISTIC_side
        if side in 'NS':
            np.testing.assert_array_almost_equal(R_v[i][0], 0)
            np.testing.assert_array_almost_equal(R_v[i][1] - w1(t+t2, *xyz), 0)
            np.testing.assert_array_almost_equal(R_v[i][2] - w2(t+t2, *xyz), 0)
        elif side in 'WE':
            np.testing.assert_array_almost_equal(R_v[i][0] - w0(t+t2, *xyz), 0)
            np.testing.assert_array_almost_equal(R_v[i][1], 0)
            np.testing.assert_array_almost_equal(R_v[i][2] - w2(t+t2, *xyz), 0)
        elif side in 'BF':
            np.testing.assert_array_almost_equal(R_v[i][0] - w0(t+t2, *xyz), 0)
            np.testing.assert_array_almost_equal(R_v[i][1] - w1(t+t2, *xyz), 0)
            np.testing.assert_array_almost_equal(R_v[i][2], 0)
        else:
            raise Exception()

    #now we test norm component + parallel component = the vector on all trace elements ------------------
    T_para.current_time = t + t1
    N_xyz, N_v = N.reconstruct(xi, et, sg) # on all trace-elements
    P_xyz, P_v = T_para.reconstruct(xi, et, sg) # on all trace-elements
    for i in N_xyz:
        xyz = N_xyz[i]
        nv = N_v[i]
        pv = P_v[i]
        vec = (nv[0] + pv[0], nv[1] + pv[1], nv[2] + pv[2])
        VEC = (w0(t+t1, *xyz), w1(t+t1, *xyz), w2(t+t1, *xyz))
        np.testing.assert_array_almost_equal(vec[0]-VEC[0], 0)
        np.testing.assert_array_almost_equal(vec[1]-VEC[1], 0)
        np.testing.assert_array_almost_equal(vec[2]-VEC[2], 0)

    # T_perp component -------------------------------------------------------------------------
    T_perp = W.components.T_perp
    T_perp.current_time = t + t3
    assert T_perp.current_time == W.current_time + t3, f"T_perp and W has decoupled!"
    assert T_perp.standard_properties.name == "T-perp-component-of-" + W.standard_properties.name, f"The naming rule."
    R_xyz, R_v = T_perp.reconstruct(xi, et, sg, i='on_mesh_boundaries')
    for i in R_xyz:
        xyz = R_xyz[i]
        te = T_perp.mesh.trace.elements[i]
        side = te.CHARACTERISTIC_side
        if side in 'NS':
            np.testing.assert_array_almost_equal(R_v[i][0], 0)
            np.testing.assert_array_almost_equal(R_v[i][1] - w2(t+t3, *xyz), 0)
            np.testing.assert_array_almost_equal(R_v[i][2] + w1(t+t3, *xyz), 0)
        elif side in 'WE':
            np.testing.assert_array_almost_equal(R_v[i][0] + w2(t+t3, *xyz), 0)
            np.testing.assert_array_almost_equal(R_v[i][1], 0)
            np.testing.assert_array_almost_equal(R_v[i][2] - w0(t+t3, *xyz), 0)
        elif side in 'BF':
            np.testing.assert_array_almost_equal(R_v[i][0] - w1(t+t3, *xyz), 0)
            np.testing.assert_array_almost_equal(R_v[i][1] + w0(t+t3, *xyz), 0)
            np.testing.assert_array_almost_equal(R_v[i][2], 0)
        else:
            raise Exception()

    # ---------------- generate a new random mesh ----------------------------------------------
    FC = random_3D_FormCaller_of_total_load_around(load, exclude_periodic=True)
    W = FC('vector', (w0, w1, w2))
    U = FC('vector', (u0, u1, u2))

    # ---------- neg ---------------------------------------------------------------------------
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

    # ----------- sub ---------------------------------------------------
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

    # ----------- add ---------------------------------------------------
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

    # ------ cross product ----------------------------------------------
    X = W.do.cross_product(U)
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
    xi = np.linspace(-0.9-random.random()/10, 0.9+random.random()/10, I)
    et = np.linspace(-0.8-random.random()/5, 0.8+random.random()/5, J)
    sg = np.linspace(-0.9-random.random()/10, 0.9+random.random()/10, K)
    x, y, z = np.meshgrid(xi, et, sg, indexing='ij')

    if rAnk == mAster_rank:
        load = random.randint(100, 499)
    else:
        load= None
    load = cOmm.bcast(load, root=mAster_rank)
    FC = random_3D_FormCaller_of_total_load_around(load, exclude_periodic=False)

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



def test_Form_NO3_3dCSCG_TensorField():
    """"""
    if rAnk == mAster_rank:
        print(f"-T- [test_Form_NO3_3dCSCG_TensorField]...", flush=True)

    t = random.random() * 10
    I, J, K = random.randint(2,10), random.randint(4,8), random.randint(3,9)
    xi = np.linspace(-0.9-random.random()/10, 0.9+random.random()/10, I)
    et = np.linspace(-0.8-random.random()/5, 0.8+random.random()/5, J)
    sg = np.linspace(-0.9-random.random()/10, 0.9+random.random()/10, K)
    x, y, z = np.meshgrid(xi, et, sg, indexing='ij')

    def T00(t, x, y, z): return 2 * t * np.sin(np.pi * x) * np.cos(3 * np.pi * y) * np.cos(np.pi * z)
    def T01(t, x, y, z): return 1.2 * t * np.sin(2 * np.pi * x) * np.cos(np.pi * y) * np.sin(np.pi * z)
    def T02(t, x, y, z): return 1.2 * t * np.sin(2 * np.pi * x) * np.sin(np.pi * y) * np.cos(np.pi * z)

    def T10(t, x, y, z): return 2 * t * np.cos(np.pi * x) * np.cos(0.11 * np.pi * y) * np.cos(np.pi * z)
    def T11(t, x, y, z): return t * np.cos(2 * np.pi * x) * np.cos(0.5 * np.pi * y) * np.sin(3 * np.pi * z)
    def T12(t, x, y, z): return t * np.cos(2 * np.pi * x) * np.sin(np.pi * y) * np.cos(2 * np.pi * z)

    def T20(t, x, y, z): return t * np.sin(np.pi * x) * np.cos(np.pi * y) * np.sin(np.pi * z)
    def T21(t, x, y, z): return t * np.sin(2 * np.pi * x) * np.cos(np.pi * y) * np.sin(0.1 * np.pi * z)
    def T22(t, x, y, z): return t * np.sin(2 * np.pi * x) * np.sin(np.pi * y) * np.sin(np.pi * z)

    def t00(t, x, y, z): return 2 * t * np.sin(np.pi * x) * np.cos(0.1 * np.pi * y) * np.cos(np.pi * z)
    def t01(t, x, y, z): return t * np.sin(5 * np.pi * x) * np.cos(0.1 * np.pi * y) * np.sin(np.pi * z)
    def t02(t, x, y, z): return t * np.sin(5 * np.pi * x) * np.sin(0.1 * np.pi * y) * np.cos(0.5 * np.pi * z)

    def t10(t, x, y, z): return 3 * t * np.cos(np.pi * x) * np.cos(np.pi * y) * np.cos(np.pi * z)
    def t11(t, x, y, z): return 1.2 * t * np.cos(3 * np.pi * x) * np.cos(np.pi * y) * np.sin(0.1 * np.pi * z)
    def t12(t, x, y, z): return t * np.cos(3 * np.pi * x) * np.sin(np.pi * y) * np.cos(np.pi * z)

    def t20(t, x, y, z): return 0.5 * t * np.sin(2*np.pi * x) * np.cos(0.1 * np.pi * y) * np.sin(np.pi * z)
    def t21(t, x, y, z): return 1.5 * t * np.sin(4 * np.pi * x) * np.cos(0.1 * np.pi * y) * np.sin(np.pi * z)
    def t22(t, x, y, z): return 0.2 * t * np.sin(4 * np.pi * x) * np.sin(np.pi * y) * np.sin(3*np.pi * z)

    if rAnk == mAster_rank:
        load = random.randint(100, 499)
    else:
        load= None
    load = cOmm.bcast(load, root=mAster_rank)
    FC = random_3D_FormCaller_of_total_load_around(load, exclude_periodic=False)

    w = FC('tensor', ([T00, T01, T02], [T10, T11, T12], [T20, T21, T22]))
    u = FC('tensor', ([t00, t01, t02], [t10, t11, t12], [t20, t21, t22]))

    #---- neg --------------------------------------------------------------------
    X = -w
    X.current_time = t
    R_xyz, R_v = X.reconstruct(x, y, z)
    for i in R_xyz:
        xyz = R_xyz[i]
        A00, A01, A02 = -T00(t, *xyz), -T01(t, *xyz), -T02(t, *xyz)
        A10, A11, A12 = -T10(t, *xyz), -T11(t, *xyz), -T12(t, *xyz)
        A20, A21, A22 = -T20(t, *xyz), -T21(t, *xyz), -T22(t, *xyz)
        C0, C1, C2 = R_v[i]
        C00, C01, C02 = C0
        C10, C11, C12 = C1
        C20, C21, C22 = C2
        assert np.max(np.abs(C00 - A00)) < 1e-10, f"neg is not accurate enough."
        assert np.max(np.abs(C01 - A01)) < 1e-10, f"neg is not accurate enough."
        assert np.max(np.abs(C02 - A02)) < 1e-10, f"neg is not accurate enough."
        assert np.max(np.abs(C10 - A10)) < 1e-10, f"neg is not accurate enough."
        assert np.max(np.abs(C11 - A11)) < 1e-10, f"neg is not accurate enough."
        assert np.max(np.abs(C12 - A12)) < 1e-10, f"neg is not accurate enough."
        assert np.max(np.abs(C20 - A20)) < 1e-10, f"neg is not accurate enough."
        assert np.max(np.abs(C21 - A21)) < 1e-10, f"neg is not accurate enough."
        assert np.max(np.abs(C22 - A22)) < 1e-10, f"neg is not accurate enough."

    #---- sub --------------------------------------------------------------------
    X = w - u
    X.current_time = t
    R_xyz, R_v = X.reconstruct(x, y, z)
    for i in R_xyz:
        xyz = R_xyz[i]
        A00, A01, A02 = T00(t, *xyz)-t00(t, *xyz), T01(t, *xyz)-t01(t, *xyz), T02(t, *xyz)-t02(t, *xyz)
        A10, A11, A12 = T10(t, *xyz)-t10(t, *xyz), T11(t, *xyz)-t11(t, *xyz), T12(t, *xyz)-t12(t, *xyz)
        A20, A21, A22 = T20(t, *xyz)-t20(t, *xyz), T21(t, *xyz)-t21(t, *xyz), T22(t, *xyz)-t22(t, *xyz)
        C0, C1, C2 = R_v[i]
        C00, C01, C02 = C0
        C10, C11, C12 = C1
        C20, C21, C22 = C2
        assert np.max(np.abs(C00 - A00)) < 1e-10, f"sub is not accurate enough."
        assert np.max(np.abs(C01 - A01)) < 1e-10, f"sub is not accurate enough."
        assert np.max(np.abs(C02 - A02)) < 1e-10, f"sub is not accurate enough."
        assert np.max(np.abs(C10 - A10)) < 1e-10, f"sub is not accurate enough."
        assert np.max(np.abs(C11 - A11)) < 1e-10, f"sub is not accurate enough."
        assert np.max(np.abs(C12 - A12)) < 1e-10, f"sub is not accurate enough."
        assert np.max(np.abs(C20 - A20)) < 1e-10, f"sub is not accurate enough."
        assert np.max(np.abs(C21 - A21)) < 1e-10, f"sub is not accurate enough."
        assert np.max(np.abs(C22 - A22)) < 1e-10, f"sub is not accurate enough."

    #---- add --------------------------------------------------------------------
    X = w + u
    X.current_time = t
    R_xyz, R_v = X.reconstruct(x, y, z)
    for i in R_xyz:
        xyz = R_xyz[i]
        A00, A01, A02 = T00(t, *xyz)+t00(t, *xyz), T01(t, *xyz)+t01(t, *xyz), T02(t, *xyz)+t02(t, *xyz)
        A10, A11, A12 = T10(t, *xyz)+t10(t, *xyz), T11(t, *xyz)+t11(t, *xyz), T12(t, *xyz)+t12(t, *xyz)
        A20, A21, A22 = T20(t, *xyz)+t20(t, *xyz), T21(t, *xyz)+t21(t, *xyz), T22(t, *xyz)+t22(t, *xyz)
        C0, C1, C2 = R_v[i]
        C00, C01, C02 = C0
        C10, C11, C12 = C1
        C20, C21, C22 = C2
        assert np.max(np.abs(C00 - A00)) < 1e-10, f"add is not accurate enough."
        assert np.max(np.abs(C01 - A01)) < 1e-10, f"add is not accurate enough."
        assert np.max(np.abs(C02 - A02)) < 1e-10, f"add is not accurate enough."
        assert np.max(np.abs(C10 - A10)) < 1e-10, f"add is not accurate enough."
        assert np.max(np.abs(C11 - A11)) < 1e-10, f"add is not accurate enough."
        assert np.max(np.abs(C12 - A12)) < 1e-10, f"add is not accurate enough."
        assert np.max(np.abs(C20 - A20)) < 1e-10, f"add is not accurate enough."
        assert np.max(np.abs(C21 - A21)) < 1e-10, f"add is not accurate enough."
        assert np.max(np.abs(C22 - A22)) < 1e-10, f"add is not accurate enough."

    return 1







if __name__ == '__main__':
    # mpiexec -n 6 python _3dCSCG\tests\unittests\fields.py
    test_Form_NO1_3dCSCG_VectorField()