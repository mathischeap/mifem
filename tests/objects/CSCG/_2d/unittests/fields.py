# -*- coding: utf-8 -*-

import sys
if './' not in sys.path:
    sys.path.append('./')
from root.config.main import *
from __init__ import cscg2
import random
from tests.objects.CSCG._2d.randObj.form_caller import random_FormCaller_of_total_load_around
from tests.objects.CSCG._2d.randObj.field import random_scalar, random_vector


def test_Fields_NO1_vector():
    """"""
    if RANK == MASTER_RANK:
        load = random.randint(100, 1000)
        print(f"~~~ [test_Fields_NO1_vector] @ load= {load}... ", flush=True)
    else:
        load = None

    load = COMM.bcast(load, root=MASTER_RANK)
    FC = random_FormCaller_of_total_load_around(load)

    I, J = random.randint(2, 10), random.randint(4, 8)
    x = np.linspace(-0.9-random.random()/10, 0.9+random.random()/10, I)
    y = np.linspace(-0.8-random.random()/5, 0.8+random.random()/5, J)
    x, y = np.meshgrid(x, y, indexing='ij')
    t = random.random()

    def u(t, x, y): return - np.pi * np.sin(2.56*np.pi*x) * np.cos(3.12*np.pi*y) * np.sin(t)/1.554
    def v(t, x, y): return np.pi * np.cos(2.112*np.pi*x) * np.sin(1.98*np.pi*y) * np.sin(t)*1.23
    def du_dt(t, x, y): return - np.pi * np.sin(2.56*np.pi*x) * np.cos(3.12*np.pi*y) * np.cos(t)/1.554
    def dv_dt(t, x, y): return np.pi * np.cos(2.112*np.pi*x) * np.sin(1.98*np.pi*y) * np.cos(t)*1.23

    v = FC('vector', (u, v))
    # test time derivatives ---------------------------------------------------
    v_t = v.numerical.time_derivative
    v_t.current_time = t
    R_xyz, R_v = v_t.reconstruct(x, y)
    for i in R_xyz:
        assert np.max(np.abs(R_v[i][0] - du_dt(t, *R_xyz[i]))) < 1e-7
        assert np.max(np.abs(R_v[i][1] - dv_dt(t, *R_xyz[i]))) < 1e-7

    #  ---- test numerical curl -------------------------------------------------------------
    def dv_dx__m__du_dy(t, x, y):
        return -2.112*np.pi**2 * np.sin(2.112*np.pi*x) * np.sin(1.98*np.pi*y) * np.sin(t)*1.23 - \
               3.12 * np.pi**2 * np.sin(2.56 * np.pi * x) * np.sin(3.12 * np.pi * y) * np.sin(t)/1.554
    rot_v = v.numerical.rot
    rot_v.current_time = t
    R_xyz, R_v = rot_v.reconstruct(x, y)
    for i in R_xyz:
        assert np.max(np.abs(R_v[i][0] - dv_dx__m__du_dy(t, *R_xyz[i]))) < 1e-7

    #  ---- test numerical div -------------------------------------------------------------
    def du_dx__p__dv_dy(t, x, y):
        return - 2.56*np.pi**2 * np.cos(2.56*np.pi*x) * np.cos(3.12*np.pi*y) * np.sin(t)/1.554 + \
               1.98*np.pi**2 * np.cos(2.112*np.pi*x) * np.cos(1.98*np.pi*y) * np.sin(t)*1.23
    div_v = v.numerical.div
    div_v.current_time = t
    R_xyz, R_v = div_v.reconstruct(x, y)
    for i in R_xyz:
        assert np.max(np.abs(R_v[i][0] - du_dx__p__dv_dy(t, *R_xyz[i]))) < 1e-7

    return 1


def test_Fields_NO2_scalar():
    """"""
    if RANK == MASTER_RANK:
        load = random.randint(100, 1000)
        print(f"~~~ [test_Fields_NO2_scalar] @ load= {load}... ", flush=True)
    else:
        load = None

    load = COMM.bcast(load, root=MASTER_RANK)
    FC = random_FormCaller_of_total_load_around(load)

    I, J = random.randint(2, 10), random.randint(4, 8)
    x = np.linspace(-0.9-random.random()/10, 0.9+random.random()/10, I)
    y = np.linspace(-0.8-random.random()/5, 0.8+random.random()/5, J)
    x, y = np.meshgrid(x, y, indexing='ij')
    t = random.random()
    def f(t, x, y): return - np.pi * np.sin(2.56*np.pi*x) * np.cos(3.12*np.pi*y) * np.sin(t) / 1.554
    def df_dx(t, x, y): return - 2.56 * np.pi**2 * np.cos(2.56*np.pi*x) * np.cos(3.12*np.pi*y) * np.sin(t) / 1.554
    def df_dy(t, x, y): return 3.12 * np.pi**2 * np.sin(2.56*np.pi*x) * np.sin(3.12*np.pi*y) * np.sin(t) / 1.554
    def df_dt(t, x, y): return - np.pi * np.sin(2.56*np.pi*x) * np.cos(3.12*np.pi*y) * np.cos(t) / 1.554
    w = FC('scalar', f)

    # test time derivatives ---------------------------------------------------
    w_t = w.numerical.time_derivative
    w_t.current_time = t
    R_xyz, R_v = w_t.reconstruct(x, y)
    for i in R_xyz:
        assert np.max(np.abs(R_v[i][0] - df_dt(t, *R_xyz[i]))) < 1e-7

    # ---------- test numerical gradient ------------------------------------------
    grad_w = w.numerical.grad
    grad_w.current_time = t
    R_xyz, R_v = grad_w.reconstruct(x, y)
    for i in R_xyz:
        assert np.max(np.abs(R_v[i][0] - df_dx(t, *R_xyz[i]))) < 1e-7
        assert np.max(np.abs(R_v[i][1] - df_dy(t, *R_xyz[i]))) < 1e-7

    # ---------- test numerical curl ------------------------------------------
    curl_w = w.numerical.curl
    curl_w.current_time = t
    R_xyz, R_v = curl_w.reconstruct(x, y)
    for i in R_xyz:
        assert np.max(np.abs(R_v[i][0] - df_dy(t, *R_xyz[i]))) < 1e-7
        assert np.max(np.abs(R_v[i][1] + df_dx(t, *R_xyz[i]))) < 1e-7

    return 1


def test_identities():
    """"""
    if RANK == MASTER_RANK:
        load = random.randint(100, 1000)
        print(f"-I- [test_identities]...", flush=True)
    else:
        load = None

    load = COMM.bcast(load, root=MASTER_RANK)
    FC = random_FormCaller_of_total_load_around(load)

    t = 1 + random.random() * 10
    I, J = random.randint(2, 10), random.randint(4, 8)
    xi = np.linspace(-0.9-random.random()/10, 0.9+random.random()/10, I)
    et = np.linspace(-0.8-random.random()/5, 0.8+random.random()/5, J)
    x, y = np.meshgrid(xi, et, indexing='ij')

    def func(t, x, y): return t * np.cos(2 * np.pi * x) * np.sin(2*np.pi * y)
    def U(t, x, y): return - 2 * t * np.sin(2*np.pi * x) * np.cos(2*np.pi * y)
    def V(t, x, y): return - t * np.cos(2 * np.pi * x) * np.sin(2*np.pi * y)

    Sca = FC('scalar', func)
    Vec = FC('vector', [U, V])
    # -------------------- test a float * a scalar or a float * a vector ----------------------------
    f = 2 + random.random()
    fS = f * Sca
    Sf = Sca * f
    fV = f * Vec
    Vf = Vec * f

    fS.current_time = t
    Sf.current_time = t
    fV.current_time = t
    Vf.current_time = t
    R_xy, R_v = fS.reconstruct(x, y)
    ____, R_V = Sf.reconstruct(x, y)
    for i in R_xy:
        xy = R_xy[i]
        f_S = f * func(t, *xy)
        np.testing.assert_array_almost_equal(R_v[i][0], f_S)
        np.testing.assert_array_almost_equal(R_V[i][0], f_S)
    R_xy, R_v = fV.reconstruct(x, y)
    ____, R_V = Vf.reconstruct(x, y)
    for i in R_xy:
        xy = R_xy[i]
        f_U = f * U(t, *xy)
        f_V = f * V(t, *xy)
        np.testing.assert_array_almost_equal(R_v[i][0], f_U)
        np.testing.assert_array_almost_equal(R_v[i][1], f_V)
        np.testing.assert_array_almost_equal(R_V[i][0], f_U)
        np.testing.assert_array_almost_equal(R_V[i][1], f_V)

    # -------------------------------- test add, sub and neg ---------------------------------------

    rS0 = random_scalar(FC.mesh)
    rS1 = random_scalar(FC.mesh)

    f0 = rS0.func[0]
    f1 = rS1.func[0]

    addS = rS0 + rS1
    subS = rS0 - rS1
    negS = - rS0

    t = random.random()
    addS.current_time = t
    R_xy, R_v = addS.reconstruct(x, y)
    for i in R_xy:
        xy = R_xy[i]
        true = f0(t, *xy) + f1(t, *xy)
        np.testing.assert_array_equal(R_v[i][0], true)
    subS.current_time = t
    R_xy, R_v = subS.reconstruct(x, y)
    for i in R_xy:
        xy = R_xy[i]
        true = f0(t, *xy) - f1(t, *xy)
        np.testing.assert_array_equal(R_v[i][0], true)
    negS.current_time = t
    R_xy, R_v = negS.reconstruct(x, y)
    for i in R_xy:
        xy = R_xy[i]
        true = -f0(t, *xy)
        np.testing.assert_array_equal(R_v[i][0], true)

    rV0 = random_vector(FC.mesh)
    rV1 = random_vector(FC.mesh)
    f0, f1 = rV0.func
    f2, f3 = rV1.func

    addV = rV0 + rV1
    subV = rV0 - rV1
    negV = - rV0
    addV.current_time = t
    R_xy, R_v = addV.reconstruct(x, y)
    for i in R_xy:
        xy = R_xy[i]
        true0 = f0(t, *xy) + f2(t, *xy)
        true1 = f1(t, *xy) + f3(t, *xy)
        np.testing.assert_array_equal(R_v[i][0], true0)
        np.testing.assert_array_equal(R_v[i][1], true1)
    subV.current_time = t
    R_xy, R_v = subV.reconstruct(x, y)
    for i in R_xy:
        xy = R_xy[i]
        true0 = f0(t, *xy) - f2(t, *xy)
        true1 = f1(t, *xy) - f3(t, *xy)
        np.testing.assert_array_equal(R_v[i][0], true0)
        np.testing.assert_array_equal(R_v[i][1], true1)
    negV.current_time = t
    R_xy, R_v = negV.reconstruct(x, y)
    for i in R_xy:
        xy = R_xy[i]
        true0 = - f0(t, *xy)
        true1 = - f1(t, *xy)
        np.testing.assert_array_equal(R_v[i][0], true0)
        np.testing.assert_array_equal(R_v[i][1], true1)

    # ---- test a scalar * a vector -------------------------------------------------------------------
    SV = Sca * Vec
    VS = Vec * Sca
    SV.current_time = t
    VS.current_time = t
    R_xy, R_v = SV.reconstruct(x, y)
    ____, R_V = VS.reconstruct(x, y)
    for i in R_xy:
        xy = R_xy[i]
        f_U = func(t, *xy) * U(t, *xy)
        f_V = func(t, *xy) * V(t, *xy)
        np.testing.assert_array_almost_equal(R_v[i][0], f_U)
        np.testing.assert_array_almost_equal(R_v[i][1], f_V)
        np.testing.assert_array_almost_equal(R_V[i][0], f_U)
        np.testing.assert_array_almost_equal(R_V[i][1], f_V)

    # -----------------------------------------------------------------------------------------------
    mesh = cscg2.mesh('crazy_periodic', c=0.)(element_layout=[5, 5])
    space = cscg2.space('polynomials')([3, 3])
    FC = cscg2.form(mesh, space)

    def W(t, x, y): return t * np.cos(2 * np.pi * x) * np.sin(2*np.pi * y)
    def U(t, x, y): return - 2 * t * np.sin(2*np.pi * x) * np.cos(2*np.pi * y)
    def V(t, x, y): return - t * np.cos(2 * np.pi * x) * np.sin(2*np.pi * y)
    def Q(t, x, y): return t * np.sin(2*np.pi * x) * np.sin(2*np.pi * y)

    w = FC('scalar', W)
    u = FC('vector', [U, V])
    q = FC('scalar', Q)

    wXu = w.do.cross_product(u)

    curl_q = q.numerical.curl

    LEFT = wXu.do.inner_product(curl_q)
    LEFT.current_time = 1
    LEFT = LEFT.do.compute_Ln_norm(n=1)

    R1 = w.do.inner_product((q * u).numerical.div)
    R2 = q.do.inner_product((w * u).numerical.div)

    R1.current_time = 1
    R2.current_time = 1
    R1 = 0.5 * R1.do.compute_Ln_norm(n=1)
    R2 = 0.5 * R2.do.compute_Ln_norm(n=1)

    np.testing.assert_almost_equal(LEFT, -R1 + R2)

    return 1


if __name__ == '__main__':
    # mpiexec -n 4 python objects/CSCG/_2d/__tests__/unittests/fields.py
    test_identities()
