# -*- coding: utf-8 -*-
import sys
if './' not in sys.path:
    sys.path.append('./')
from root.config.main import *
from tests.objects.CSCG._2d.randObj.form_caller import random_FormCaller_of_total_load_around
from objects.CSCG._2d.master import MeshGenerator, SpaceInvoker, FormCaller
from scipy.sparse import linalg as spspalinalg
import random
from tools.miLinearAlgebra.linearSystem.main import LinearSystem
from tools.elementwiseCache.dataStructures.operators.bmat.main import bmat
from tools.elementwiseCache.dataStructures.operators.concatenate.main import concatenate


def test_Form_NO1_coboundary():
    """"""
    if RANK == MASTER_RANK:
        print("~~~ [test_Form_NO1_coboundary] ...... ", flush=True)

    mesh = MeshGenerator('crazy_periodic', c=0.1)([8, 9])
    space = SpaceInvoker('polynomials')([('Lobatto', 6), ('Lobatto', 5)])
    FC = FormCaller(mesh, space)

    def p(t, x, y): return np.cos(2*np.pi*x) * np.cos(2*np.pi*y) + t/2

    def p_dpdx_i(t, x, y): return -2*np.pi * np.sin(2*np.pi*x) * np.cos(2*np.pi*y) + 0 * t
    def p_dpdy_i(t, x, y): return -2*np.pi * np.cos(2*np.pi*x) * np.sin(2*np.pi*y) + 0 * t

    def P(t, x, y): return - np.pi * np.sin(2*np.pi*x) * np.cos(2*np.pi*y) + t/2
    def Q(t, x, y): return np.pi * np.cos(2*np.pi*x) * np.sin(2*np.pi*y) + t

    def dQdx_m_dPdy(t, x, y):   # curl of a 2-d vector
        dQdx = - 2 * np.pi**2 * np.sin(2*np.pi*x) * np.sin(2*np.pi*y) + 0 * t
        dPdy = 2 * np.pi**2 * np.sin(2*np.pi*x) * np.sin(2*np.pi*y) + 0 * t
        return dQdx - dPdy

    def p_dpdy(t, x, y): return -2*np.pi * np.cos(2*np.pi*x) * np.sin(2*np.pi*y) + 0 * t
    def m_dpdx(t, x, y): return 2*np.pi * np.sin(2*np.pi*x) * np.cos(2*np.pi*y) + 0 * t

    def Po(t, x, y): return np.pi * np.sin(2*np.pi*x) * np.cos(2*np.pi*y) + t/2
    def Qo(t, x, y): return np.pi * np.cos(2*np.pi*x) * np.sin(2*np.pi*y) + t

    def DIV_PQ(t, x, y):
        dPdx = 2*np.pi**2 * np.cos(2*np.pi*x) * np.cos(2*np.pi*y) + 0 * t
        dQdy = 2*np.pi**2 * np.cos(2*np.pi*x) * np.cos(2*np.pi*y) + 0 * t
        return dPdx + dQdy
    
    sp = FC('scalar', p)
    if RANK == MASTER_RANK:
        T = [random.uniform(-10, 10) for _ in range(3)]
    else:
        T = None

    T = COMM.bcast(T, root=MASTER_RANK)

    # inner outer ------------------------------------------------
    for t in T:
        # test inner -------------------------------------------
        f0 = FC('0-f-i', hybrid=False)
        f0.CF = sp
        f0.CF.current_time = t
        f0.discretize()
        assert f0.error.L(upon=False) < 9e-7

        f1 = f0.coboundary()

        vector = FC('vector', (p_dpdx_i, p_dpdy_i))
        f1.CF = vector
        f1.CF.current_time = t
        assert f1.error.L(upon=False) < 1e-4

        f2 = f1.coboundary()
        zero = FC('scalar', 0)
        f2.CF = zero
        f2.CF.current_time = t
        np.testing.assert_almost_equal(f2.error.L(upon=False), 0)

        vector = FC('vector', (P, Q))
        f1.CF = vector
        f1.CF.current_time = t
        f1.discretize()
        assert f1.error.L(upon=False) < 2e-5

        f2 = f1.coboundary()
        scalar = FC('scalar', dQdx_m_dPdy)
        f2.CF = scalar
        f2.CF.current_time = t
        assert f2.error.L(upon=False) < 5e-4

        # test outer -------------------------------------------

        f0 = FC('0-f-o', hybrid=False)
        f0.CF = sp
        f0.CF.current_time = t
        f0.discretize()
        assert f0.error.L(upon=False) < 9e-7

        f1 = f0.coboundary()

        vector = FC('vector', (p_dpdy, m_dpdx))
        f1.CF = vector
        f1.CF.current_time = t

        assert f1.error.L(upon=False) < 1e-4

        f2 = f1.coboundary()
        zero = FC('scalar', 0)
        f2.CF = zero
        f2.CF.current_time = t
        np.testing.assert_almost_equal(f2.error.L(upon=False), 0)

        vector = FC('vector', (Po, Qo))
        f1.CF = vector
        f1.CF.current_time = t
        f1.discretize()
        assert f1.error.L(upon=False) < 3e-5

        f2 = f1.coboundary()
        scalar = FC('scalar', DIV_PQ)
        f2.CF = scalar
        f2.CF.current_time = t
        assert f2.error.L(upon=False) < 5e-4

    return 1


def test_Form_NO2_Naive_numbering():
    """"""
    if RANK == MASTER_RANK:
        print("~~~ [test_Form_NO2_Naive_numbering] ...... ", flush=True)

    mesh = MeshGenerator('crazy', c=0.1)([3, 2], EDM='debug')
    space = SpaceInvoker('polynomials')([('Lobatto', 2), ('Lobatto', 3)])
    FC = FormCaller(mesh, space)

    f0 = FC('0-f-i', hybrid=True)
    f1i = FC('1-f-i', hybrid=True)
    f1o = FC('1-f-o', hybrid=True)
    f2 = FC('2-f-i', hybrid=True)

    t1o = FC('1-t-o')

    G0 = f0.numbering.gathering
    G1I = f1i.numbering.gathering
    G1O = f1o.numbering.gathering
    G2 = f2.numbering.gathering
    T1O = t1o.numbering.gathering

    for i in mesh.elements:
        np.testing.assert_array_equal(
            G0[i].full_vector, np.arange(i * f0.num.basis, (i+1) * f0.num.basis))
        np.testing.assert_array_equal(
            G1I[i].full_vector, np.arange(i * f1i.num.basis, (i+1) * f1i.num.basis))
        np.testing.assert_array_equal(
            G1O[i].full_vector, np.arange(i * f1o.num.basis, (i+1) * f1o.num.basis))
        np.testing.assert_array_equal(
            G2[i].full_vector, np.arange(i * f2.num.basis, (i+1) * f2.num.basis))

    if 1 in mesh.elements:
        np.testing.assert_array_equal(
            T1O[1].full_vector,
            np.array([3, 4, 5, 10, 11, 12, 13, 14, 15, 16])
        )
    if 4 in mesh.elements:
        np.testing.assert_array_equal(
            T1O[4].full_vector,
            np.array([27, 28, 29, 32, 33, 34, 15, 16, 35, 36])
        )

    f0 = FC('0-f-i', hybrid=False)
    f1i = FC('1-f-i', hybrid=False)
    f1o = FC('1-f-o', hybrid=False)
    f2 = FC('2-f-i', hybrid=False)

    G0 = f0.numbering.gathering
    G1I = f1i.numbering.gathering
    G1O = f1o.numbering.gathering
    G2 = f2.numbering.gathering

    if 1 in mesh.elements:
        np.testing.assert_array_equal(
            G0[1].full_vector,
            np.array([2, 12, 16, 5, 13, 17, 8, 14, 18, 11, 15, 19])
        )
    if 5 in mesh.elements:
        np.testing.assert_array_equal(
            G0[5].full_vector,
            np.array([19, 23, 27, 40, 43, 46, 41, 44, 47, 42, 45, 48])
        )

    if 1 in mesh.elements:
        np.testing.assert_array_equal(
            G1I[1].full_vector,
            np.array([17, 21, 18, 22, 19, 23, 20, 24, 10, 25, 28, 13, 26, 29, 16, 27, 30])
        )
    if 3 in mesh.elements:
        np.testing.assert_array_equal(
            G1I[3].full_vector,
            np.array([6, 7, 45, 48, 46, 49, 47, 50, 51, 54, 57, 52, 55, 58, 53, 56, 59])
        )

    if 2 in mesh.elements:
        np.testing.assert_array_equal(
            G1O[2].full_vector,
            np.array([20, 31, 34, 21, 32, 35, 22, 33, 36, 37, 41, 38, 42, 39, 43, 40, 44])
        )
    if 4 in mesh.elements:
        np.testing.assert_array_equal(
            G1O[4].full_vector,
            np.array([51, 60, 63, 52, 61, 64, 53, 62, 65, 26, 30, 66, 69, 67, 70, 68, 71])
        )

    for i in mesh.elements:
        np.testing.assert_array_equal(
            G2[i].full_vector, np.arange(i * f2.num.basis, (i+1) * f2.num.basis))

    return 1


def test_Form_NO3_mass_matrices():
    """"""
    if RANK == MASTER_RANK:
        print("~~~ [test_Form_NO3_mass_matrices] ...... ", flush=True)

    mesh = MeshGenerator('crazy', c=0.12, bounds=((0, 2), (0, 3)))([2, 3], EDM='debug')
    space = SpaceInvoker('polynomials')([('Lobatto', 3), ('Lobatto', 2)])
    FC = FormCaller(mesh, space)

    f0 = FC('0-f-i', hybrid=False)
    f1i = FC('1-f-i', hybrid=False)
    f1o = FC('1-f-o', hybrid=False)
    f2 = FC('2-f-i', hybrid=False)

    M0 = f0.matrices.mass
    M1I = f1i.matrices.mass
    M1O = f1o.matrices.mass
    M2 = f2.matrices.mass

    if 0 in M1I:
        np.testing.assert_array_almost_equal(
            M1I[0].toarray()[0],
            np.array([
                0.60122847, -0.07187971,  0.01084695,  0.31006408, -0.03502005,
                0.00256356, -0.12646503,  0.01871769, -0.00541012, -0.04149325,
                -0.07118773,  0.01640982, -0.00061779,  0.01947493,  0.0279004,
                -0.00217755, -0.00156191
            ])
        )

        np.testing.assert_array_almost_equal(
            M1I[0].toarray()[5],
            np.array([
                2.56355608e-03, -3.49248356e-02,  3.39292002e-01,  8.79317872e-02,
                -3.43667975e-01,  3.31907618e+00,  1.90769236e-02, -4.86792311e-02,
                4.53067783e-01, -5.71518131e-04,  1.25425544e-03,  1.22704151e-01,
                1.25499757e-01,  6.05666173e-03, -3.39883326e-02,  2.98394338e-01,
                2.64620426e-01
            ])
        )
        np.testing.assert_array_almost_equal(
            M0[0].toarray()[4],
            np.array([
                0.00554578,  0.00184104, -0.00249636,  0.00098171,  0.05000453,
                0.02104327, -0.0189073,  0.006915,  0.00648875,  0.00315156,
                -0.00204593,  0.00072488
            ])
        )
        np.testing.assert_array_almost_equal(
            M2[0].toarray()[4],
            np.array([0.1989124, -1.01161361,  0.17135838, -1.37435142,  7.3285936, -1.53333516])
        )

    if 5 in M1I:
        np.testing.assert_array_almost_equal(
            M1I[5].toarray()[16],
            np.array([
                -0.00246943,  0.00999065, -0.04423281,  0.01011988, -0.04317603,
                0.15730403,  0.00479624, -0.02344232,  0.05624134, -0.00448805,
                0.01015547, -0.00717787, -0.02391952,  0.03549893, -0.08172606,
                0.04729042,  0.17342575
            ])
        )
        np.testing.assert_array_almost_equal(
            M1I[5].toarray()[11],
            np.array([
                0.00265893,  0.03943992,  0.20932935, -0.00753896,  0.09378854,
                0.3155266,  0.00499139, -0.03020874, -0.0489097, -0.07634322,
                0.114424,  0.80176059,  0.06348914,  0.01112712, -0.02682325,
                -0.13038107, -0.00717787
            ])
        )
        np.testing.assert_array_almost_equal(
            M1O[5].toarray()[11],
            np.array([
                0.126342,  0.157021,  0.007539, -0.015498,  0.067897,  0.052683,
                0.028214, -0.01012,  0.223465, -0.039418,  0.014051,  1.739309,
                -0.30496,  0.129221,  0.233822, -0.036199,  0.01639
            ])
        )

        np.testing.assert_array_almost_equal(
            M0[5].toarray()[11],
            np.array([
                -0.00036701,  0.00063908, -0.00052651, -0.00175327,  0.0006056,
                -0.00105295,  0.00170828,  0.00397803,  0.0013024, -0.00237512,
                0.00323713,  0.00803215
            ])
        )
        np.testing.assert_array_almost_equal(
            M0[5].toarray()[7],
            np.array([
                0.00086242, -0.00150339,  0.00039775,  0.00303505,  0.00578341,
                -0.00948722,  0.00735124,  0.02618595,  0.0006056, -0.00105295,
                0.00170828,  0.00397803
            ])
        )

    if 3 in M1I:
        np.testing.assert_array_almost_equal(
            M2[3].toarray()[0],
            np.array([11.5256843, -1.47127256,  0.57304065, -1.58437845,  0.2166255, -0.08919671])
        )
        np.testing.assert_array_almost_equal(
            M2[3].toarray()[1],
            np.array([-1.47127256,  5.64509819, -1.25453527,  0.2166255, -0.9143819, 0.2166255])
        )
        np.testing.assert_array_almost_equal(
            M2[3].toarray()[2],
            np.array([0.57304065, -1.25453527,  8.15018015, -0.08919671,  0.2166255, -1.58437845])
        )
        np.testing.assert_array_almost_equal(
            M2[3].toarray()[3],
            np.array([-1.58437845,  0.2166255, -0.08919671,  8.15018015, -1.25453527, 0.57304065])
        )
        np.testing.assert_array_almost_equal(
            M2[3].toarray()[4],
            np.array([0.2166255, -0.9143819,  0.2166255, -1.25453527,  5.64509819, -1.47127256])
        )
        np.testing.assert_array_almost_equal(
            M2[3].toarray()[5],
            np.array([-0.08919671,  0.2166255, -1.58437845,  0.57304065, -1.47127256, 11.5256843])
        )

    return 1


def test_Form_NO4_Hodge():
    """"""
    if RANK == MASTER_RANK:
        print("~~~ [test_Form_NO4_Hodge] ...... ", flush=True)

    def p(t, x, y): return np.cos(2*np.pi*x) * np.cos(2*np.pi*y) + t/2

    def P(t, x, y): return - np.pi * np.sin(2*np.pi*x) * np.cos(2*np.pi*y) + t/2
    def Q(t, x, y): return np.pi * np.cos(2*np.pi*x) * np.sin(2*np.pi*y) + t

    mesh = MeshGenerator('crazy', c=0.1, bounds=((0, 1), (0, 1)))([8, 9])
    space = SpaceInvoker('polynomials')([('Lobatto', 4), ('Lobatto', 3)])
    FC = FormCaller(mesh, space)

    scalar = FC('scalar', p)
    vector = FC('vector', (P, Q))
    t = random.uniform(-1, 1)

    f0 = FC('0-f-i', hybrid=False)
    f1i = FC('1-f-i', hybrid=False)
    f1o = FC('1-f-o', hybrid=False)
    f2 = FC('2-f-i', hybrid=False)

    # ------- f0 f2 -------------------------
    f0.CF = scalar
    f0.CF.current_time = t
    f0.discretize()

    M = f2.matrices.mass
    W = f0.operators.wedge(f2)
    f2_cochain = dict()
    for i in M:
        Mi = M[i].tocsc()
        Wi = W[i]
        H = spspalinalg.inv(Mi) @ Wi
        f2_cochain[i] = H @ f0.cochain.local[i]

    f2.cochain.local = f2_cochain
    f2.CF = scalar
    f2.CF.current_time = t
    f2_error = f2.error.L(upon=False)
    assert f2_error < 0.0015

    f2.discretize()
    M = f0.matrices.mass
    W = f2.operators.wedge(f0)
    f0_cochain = dict()
    for i in M:
        Mi = M[i].tocsc()
        Wi = W[i]
        H = spspalinalg.inv(Mi) @ Wi
        # noinspection PyUnresolvedReferences
        f0_cochain[i] = H @ f2.cochain.local[i]

    f0.cochain.local = f0_cochain
    f0_error = f0.error.L(upon=False)
    assert f0_error < 0.0015

    # ---------- f1 inner f1 outer ----------------
    f1i.CF = vector
    f1i.CF.current_time = t
    f1i.discretize()
    f1i_error = f1i.error.L(upon=False)
    assert f1i_error < 0.0035

    M = f1o.matrices.mass
    W = f1i.operators.wedge(f1o)
    f1o_cochain = dict()
    for i in M:
        Mi = M[i].tocsc()
        Wi = W[i]
        H = spspalinalg.inv(Mi) @ Wi
        f1o_cochain[i] = H @ f1i.cochain.local[i]

    f1o.cochain.local = f1o_cochain
    f1o.CF = vector
    f1o.CF.current_time = t
    f1o_error = f1o.error.L(upon=False)
    assert f1o_error < 0.006

    f1o.discretize()
    M = f1i.matrices.mass
    W = f1o.operators.wedge(f1i)
    f1i_cochain = dict()
    for i in M:
        Mi = M[i].tocsc()
        Wi = W[i]
        H = spspalinalg.inv(Mi) @ Wi
        # noinspection PyUnresolvedReferences
        f1i_cochain[i] = H @ f1o.cochain.local[i]

    f1i.cochain.local = f1i_cochain
    f1i_error = f1i.error.L(upon=False)
    assert f1i_error < 0.006
    return 1


def test_Form_NO5_cross_product():
    """"""
    if RANK == MASTER_RANK:
        t = random.uniform(-1, 1)
        c = random.randint(0, 3) * random.random() / 10
        IH = [True, False][random.randint(0, 1)]
        print(f"~~~ [test_Form_NO5_cross_product] @ c= %.4f... " % c, flush=True)
    else:
        t = None
        c = None
        IH = None

    t, c, IH = COMM.bcast([t, c, IH], root=MASTER_RANK)

    def w(t, x, y): return np.cos(2*np.pi*x) * np.cos(2*np.pi*y) + t/2
    def u(t, x, y): return - np.pi * np.sin(2*np.pi*x) * np.cos(2*np.pi*y) + t/2
    def v(t, x, y): return np.pi * np.cos(2*np.pi*x) * np.sin(2*np.pi*y) + t
    def a(t, x, y): return np.cos(2.08*np.pi*x) * np.sin(4.5*np.pi*y) + t/2.336
    def b(t, x, y): return np.sin(3.1*np.pi*x) * np.cos(1.587*np.pi*y) + t*3.5123

    def wva_wub(t, x, y): return - w(t, x, y) * v(t, x, y) * a(t, x, y) + w(t, x, y) * u(t, x, y) * b(t, x, y)

    mesh = MeshGenerator('crazy', c=c)([6, 5])
    space = SpaceInvoker('polynomials')([('Lobatto', 4), ('Lobatto', 3)])
    FC = FormCaller(mesh, space)

    scalar1 = FC('scalar', w)
    vector1 = FC('vector', (u, v))
    vector2 = FC('vector', (a, b))
    scalar2 = FC('scalar', wva_wub)

    w0 = FC('0-f-o', hybrid=IH)
    u1 = FC('1-f-o', hybrid=IH)
    e1 = FC('1-f-o', hybrid=IH)
    R2 = FC('2-f-o', hybrid=IH)

    w0.CF = scalar1
    w0.CF.current_time = t
    w0.discretize()
    u1.CF = vector1
    u1.CF.current_time = t
    u1.discretize()
    e1.CF = vector2
    e1.CF.current_time = t
    e1.discretize()

    R2.CF = scalar2
    R2.CF.current_time = t
    R2.discretize()

    CP = w0.special.cross_product_1f__ip_1f(u1, e1)
    for i in mesh.elements:
        cpi = CP[i].toarray()
        result = np.einsum('j, jk, k ->', u1.cochain.local[i], cpi, u1.cochain.local[i],
                           optimize='greedy')
        np.testing.assert_almost_equal(result, 0)

    return 1


def test_Form_NO6_reconstruction_matrices():
    """"""
    if RANK == MASTER_RANK:
        load = random.randint(10, 1000)
        IH = [True, False][random.randint(0, 1)]
        print(f"~~~ [test_Form_NO6_reconstruction_matrices] @ load= {load}... ", flush=True)
    else:
        load = None
        IH = None
    load, IH = COMM.bcast([load, IH], root=MASTER_RANK)

    # FC = random_FormCaller_of_total_load_around(load, mesh_pool='rectangle_periodic')
    FC = random_FormCaller_of_total_load_around(load)

    a, b, c, d = random.random(), random.random(), random.random(), random.random()
    rT = random.random()
    def P(t, x, y): return - np.sin(a*np.pi*x) * np.cos(b*np.pi*y) + t/3
    def Q(t, x, y): return np.cos(c*np.pi*x) * np.sin(d*np.pi*y) - t/3

    scalar = FC('scalar', P)
    vector = FC('vector', (P, Q))
    xi = np.linspace(-1, 1, 11)
    eta = np.linspace(-1, 1, 9)
    # ----------- outer 0-form -------------
    f0o = FC('0-f-o', hybrid=IH)
    f0i = FC('0-f-i', hybrid=IH)
    for f0 in (f0o, f0i):
        f0.CF = scalar
        f0.CF.current_time = rT
        f0.discretize()
        R = f0.reconstruct(xi, eta, ravel=True)[1]
        RR = f0.reconstruct(xi, eta, ravel=True, vectorized=True, value_only=True)[0]
        MR = f0.do.make_reconstruction_matrix_on_grid(xi, eta)
        for _, i in enumerate(MR):
            r = MR[i] @ f0.cochain.local[i]
            np.testing.assert_almost_equal(np.max(np.abs(r - R[i][0])), 0)
            np.testing.assert_almost_equal(np.max(np.abs(r - RR[_, :])), 0)

        LnEnF = f0.do.compute_Ln_energy(vectorized=False)
        LnEnT = f0.do.compute_Ln_energy(vectorized=True)
        np.testing.assert_almost_equal(LnEnF, LnEnT, decimal=5)

    # --------- outer 1-form ---------------------
    f1o = FC('1-f-o', hybrid=IH)
    f1o.CF = vector
    f1o.CF.current_time = rT
    f1o.discretize()

    R = f1o.reconstruct(xi, eta, ravel=True)[1]
    RR = f1o.reconstruct(xi, eta, ravel=True, vectorized=True, value_only=True)
    MR = f1o.do.make_reconstruction_matrix_on_grid(xi, eta)
    for _, i in enumerate(MR):
        u, v = MR[i] @ f1o.cochain.local[i]
        U, V = R[i]
        np.testing.assert_almost_equal(np.max(np.abs(U - u)), 0, decimal=5)
        np.testing.assert_almost_equal(np.max(np.abs(v - V)), 0, decimal=5)
        U, V = RR[0][_, :], RR[1][_, :]
        np.testing.assert_almost_equal(np.max(np.abs(U - u)), 0, decimal=5)
        np.testing.assert_almost_equal(np.max(np.abs(v - V)), 0, decimal=5)

    LnEnF = f1o.do.compute_Ln_energy(vectorized=False)
    LnEnT = f1o.do.compute_Ln_energy(vectorized=True)
    np.testing.assert_almost_equal(LnEnF, LnEnT, decimal=5)
    L2Energy = f1o.do.compute_L2_energy_with()
    np.testing.assert_almost_equal(LnEnT, L2Energy, decimal=1)

    # ------------ inner 1-form ------------------------------
    vector = FC('vector', (P, Q))
    f1i = FC('1-f-i', hybrid=IH)
    f1i.CF = vector
    f1i.CF.current_time = rT
    f1i.discretize()

    R = f1i.reconstruct(xi, eta, ravel=True)[1]
    RR = f1i.reconstruct(xi, eta, ravel=True, vectorized=True, value_only=True)
    MR = f1i.do.make_reconstruction_matrix_on_grid(xi, eta)
    for _, i in enumerate(MR):
        u, v = MR[i] @ f1i.cochain.local[i]
        U, V = R[i]
        np.testing.assert_almost_equal(np.max(np.abs(U - u)), 0, decimal=5)
        np.testing.assert_almost_equal(np.max(np.abs(v - V)), 0, decimal=5)
        U, V = RR[0][_, :], RR[1][_, :]
        np.testing.assert_almost_equal(np.max(np.abs(U - u)), 0, decimal=5)
        np.testing.assert_almost_equal(np.max(np.abs(v - V)), 0, decimal=5)

    LnEnF = f1i.do.compute_Ln_energy(vectorized=False)
    LnEnT = f1i.do.compute_Ln_energy(vectorized=True)
    np.testing.assert_almost_equal(LnEnF, LnEnT, decimal=5)
    L2Energy = f1i.do.compute_L2_energy_with()
    np.testing.assert_almost_equal(LnEnT, L2Energy, decimal=1)

    # ----------- 2-form -------------
    f2o = FC('2-f-o', hybrid=IH)
    f2i = FC('2-f-i', hybrid=IH)
    for f2 in (f2o, f2i):
        f2.CF = scalar
        f2.CF.current_time = rT
        f2.discretize()
        R = f2.reconstruct(xi, eta, ravel=True)[1]
        RR = f2.reconstruct(xi, eta, ravel=True, vectorized=True, value_only=True)[0]
        MR = f2.do.make_reconstruction_matrix_on_grid(xi, eta)
        for _, i in enumerate(MR):
            r = MR[i] @ f2.cochain.local[i]
            np.testing.assert_almost_equal(np.max(np.abs(r - R[i][0])), 0)
            np.testing.assert_almost_equal(np.max(np.abs(r - RR[_, :])), 0)

        LnEnF = f2.do.compute_Ln_energy(vectorized=False)
        LnEnT = f2.do.compute_Ln_energy(vectorized=True)
        np.testing.assert_almost_equal(LnEnF, LnEnT)

        L2Energy = f2.do.compute_L2_energy_with()
        np.testing.assert_almost_equal(LnEnT, L2Energy, decimal=1)

    return 1


def test_Form_NO7_weak_curl():
    """"""
    if RANK == MASTER_RANK:
        print(f"~~~ [test_Form_NO7_weak_curl]... ", flush=True)
    else:
        pass

    # --- test 1: non-hybrid; periodic boundary condition ----------------------------------
    def u(t, x, y): return - np.pi * np.sin(2*np.pi*x) * np.cos(2*np.pi*y) + t/1.554
    def v(t, x, y): return np.pi * np.cos(2*np.pi*x) * np.sin(2*np.pi*y) + t*1.23
    mesh = MeshGenerator('crazy_periodic', c=0.1)([25, 15], EDM=None)
    space = SpaceInvoker('polynomials')([('Lobatto', 3), ('Lobatto', 4)])
    FC = FormCaller(mesh, space)
    U = FC('vector', (u, v))
    rot_U = U.numerical.rot

    w0 = FC('0-f-o', hybrid=False)  # w0 = curl (u1)
    u1 = FC('1-f-o', hybrid=False)  # w0 = curl (u1)

    u1.CF = U
    u1.CF.current_time = 1
    u1.discretize()

    M0 = w0.matrices.mass

    M1 = u1.matrices.mass
    E01 = w0.matrices.incidence.T
    b = concatenate([E01 @ M1 @ u1.cochain.EWC, ])
    A = bmat(([M0, ], ))

    A.gathering_matrices = (w0, w0)
    b.gathering_matrix = w0
    LS = LinearSystem(A, b)

    result = LS.solve('GMRES')(0, restart=100, maxiter=50)

    w0.CF = rot_U
    w0.CF.current_time = 1
    result[0].do.distributed_to(w0)

    assert w0.error.L() < 0.0015, f"periodic boundary condition test fails."

    # -------------------- 0 tangent velocity boundary condition ----------------------
    def u(t, x, y): return np.cos(2.5984*np.pi*x) * np.sin(2*np.pi*y) + t
    def v(t, x, y): return np.sin(2*np.pi*x) * np.cos(3.2158*np.pi*y) + t
    mesh = MeshGenerator('crazy', c=0.1)([25, 15], EDM=None)
    space = SpaceInvoker('polynomials')([('Lobatto', 3), ('Lobatto', 4)])
    FC = FormCaller(mesh, space)
    U = FC('vector', (u, v))
    rot_U = U.numerical.rot

    w0 = FC('0-f-o', hybrid=False)  # w0 = curl (u1)
    u1 = FC('1-f-o', hybrid=False)  # w0 = curl (u1)

    u1.CF = U
    u1.CF.current_time = 0
    u1.discretize()

    M0 = w0.matrices.mass

    M1 = u1.matrices.mass
    E01 = w0.matrices.incidence.T
    b = concatenate([E01 @ M1 @ u1.cochain.EWC, ])
    A = bmat(([M0, ], ))

    A.gathering_matrices = (w0, w0)
    b.gathering_matrix = w0
    LS = LinearSystem(A, b)

    result = LS.solve('GMRES')(0, restart=100, maxiter=50)

    w0.CF = rot_U
    w0.CF.current_time = 0
    result[0].do.distributed_to(w0)

    assert w0.error.L() < 0.0006, f"0-tangent velocity boundary condition test fails."

    return 1


if __name__ == '__main__':
    # mpiexec -n 4 python tests/objects/CSCG/_2d/unittests/standard_forms/general.py
    # test_Form_NO5_cross_product()
    # test_Form_NO6_reconstruction_matrices()
    # test_Form_NO6_reconstruction_matrices()
    test_Form_NO3_mass_matrices()
