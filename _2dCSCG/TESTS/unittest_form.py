
import sys
if './' not in sys.path: sys.path.append('./')
from root.config import *
import random
from _2dCSCG.main import MeshGenerator, SpaceInvoker, FormCaller
from scipy.sparse import linalg as spspalinalg



def test_Form_NO1_coboundary():
    """"""
    if rAnk == mAster_rank:
        print("~~~ [test_Form_NO1_coboundary] ...... ", flush=True)

    mesh = MeshGenerator('crazy_periodic', c=0.1)([8, 9])
    space = SpaceInvoker('polynomials')([('Lobatto', 6), ('Lobatto', 5)])
    FC = FormCaller(mesh, space)

    def p(t, x, y): return np.cos(2*np.pi*x) * np.cos(2*np.pi*y) + t/2

    def p_dpdx_i(t, x, y): return -2*np.pi * np.sin(2*np.pi*x) * np.cos(2*np.pi*y) + 0* t
    def p_dpdy_i(t, x, y): return -2*np.pi * np.cos(2*np.pi*x) * np.sin(2*np.pi*y) + 0* t

    def P(t, x, y): return - np.pi * np.sin(2*np.pi*x) * np.cos(2*np.pi*y) + t/2
    def Q(t, x, y): return np.pi * np.cos(2*np.pi*x) * np.sin(2*np.pi*y) + t

    def dQdx_m_dPdy(t, x, y): # curl of a 2-d vector
        dQdx = - 2 * np.pi**2 * np.sin(2*np.pi*x) * np.sin(2*np.pi*y) + 0* t
        dPdy = 2 * np.pi**2 * np.sin(2*np.pi*x) * np.sin(2*np.pi*y) + 0* t
        return dQdx - dPdy

    def p_dpdy(t, x, y): return -2*np.pi * np.cos(2*np.pi*x) * np.sin(2*np.pi*y) + 0* t
    def m_dpdx(t, x, y): return  2*np.pi * np.sin(2*np.pi*x) * np.cos(2*np.pi*y) + 0* t

    def Po(t, x, y): return np.pi * np.sin(2*np.pi*x) * np.cos(2*np.pi*y) + t/2
    def Qo(t, x, y): return np.pi * np.cos(2*np.pi*x) * np.sin(2*np.pi*y) + t

    def DIV_PQ(t, x, y):
        dPdx = 2*np.pi**2 * np.cos(2*np.pi*x) * np.cos(2*np.pi*y) + 0* t
        dQdy = 2*np.pi**2 * np.cos(2*np.pi*x) * np.cos(2*np.pi*y) + 0* t
        return dPdx + dQdy
    
    sp = FC('scalar', p)

    if rAnk == mAster_rank:
        T = [random.uniform(-10, 10) for _ in range(3)]
    else:
        T = None
    T = cOmm.bcast(T, root=mAster_rank)

    # ---------- inner outer --------------------------------------------------
    for t in T:
        f0 = FC('0-f-i', is_hybrid=False)
        f0.TW.func.___DO_set_func_body_as___(sp)
        f0.TW.___DO_push_all_to_instant___(t)
        f0.discretize()
        assert f0.error.L(upon=False) < 9e-7

        f1 = f0.coboundary()

        vector = FC('vector', (p_dpdx_i, p_dpdy_i))
        f1.TW.func.___DO_set_func_body_as___(vector)
        f1.TW.___DO_push_all_to_instant___(t)
        assert f1.error.L(upon=False) < 1e-4

        f2 = f1.coboundary()
        zero = FC('scalar', 0)
        f2.TW.func.___DO_set_func_body_as___(zero)
        f2.TW.___DO_push_all_to_instant___(t)
        np.testing.assert_almost_equal(f2.error.L(upon=False), 0)

        vector = FC('vector', (P, Q))
        f1.TW.func.body = vector
        f1.TW.___DO_push_all_to_instant___(t)
        f1.discretize()
        assert f1.error.L(upon=False) < 2e-5

        f2 = f1.coboundary()
        scalar = FC('scalar', dQdx_m_dPdy)
        f2.TW.func.___DO_set_func_body_as___(scalar)
        f2.TW.___DO_push_all_to_instant___(t)
        assert f2.error.L(upon=False) < 5e-4

        # ---------- test outer --------------------------------------------------

        f0 = FC('0-f-o', is_hybrid=False)
        f0.TW.func.___DO_set_func_body_as___(sp)
        f0.TW.___DO_push_all_to_instant___(t)
        f0.discretize()
        assert f0.error.L(upon=False) < 9e-7

        f1 = f0.coboundary()

        vector = FC('vector', (p_dpdy, m_dpdx))
        f1.TW.func.body = vector
        f1.TW.___DO_push_all_to_instant___(t)

        assert f1.error.L(upon=False) < 1e-4

        f2 = f1.coboundary()
        zero = FC('scalar', 0)
        f2.TW.func.___DO_set_func_body_as___(zero)
        f2.TW.___DO_push_all_to_instant___(t)
        np.testing.assert_almost_equal(f2.error.L(upon=False), 0)

        vector = FC('vector', (Po, Qo))
        f1.TW.func.___DO_set_func_body_as___(vector)
        f1.TW.___DO_push_all_to_instant___(t)
        f1.discretize()
        assert f1.error.L(upon=False) < 3e-5

        f2 = f1.coboundary()
        scalar = FC('scalar', DIV_PQ)
        f2.TW.func.___DO_set_func_body_as___(scalar)
        f2.TW.___DO_push_all_to_instant___(t)
        assert f2.error.L(upon=False) < 5e-4

    return 1


def test_Form_NO2_Naive_numbering():
    """"""
    if rAnk == mAster_rank:
        print("~~~ [test_Form_NO2_Naive_numbering] ...... ", flush=True)

    mesh = MeshGenerator('crazy', c=0.1)([3, 2], EDM='debug')
    space = SpaceInvoker('polynomials')([('Lobatto', 2), ('Lobatto', 3)])
    FC = FormCaller(mesh, space)

    f0 = FC('0-f-i', is_hybrid=True)
    f1i = FC('1-f-i', is_hybrid=True)
    f1o = FC('1-f-o', is_hybrid=True)
    f2 = FC('2-f-i', is_hybrid=True)

    t1o = FC('1-t-o')

    G0 = f0.numbering.gathering
    G1I = f1i.numbering.gathering
    G1O = f1o.numbering.gathering
    G2 = f2.numbering.gathering
    T1O = t1o.numbering.gathering

    for i in mesh.elements:
        np.testing.assert_array_equal(
            G0[i].full_vector, np.arange( i * f0.NUM_basis, (i+1) * f0.NUM_basis))
        np.testing.assert_array_equal(
            G1I[i].full_vector, np.arange( i * f1i.NUM_basis, (i+1) * f1i.NUM_basis))
        np.testing.assert_array_equal(
            G1O[i].full_vector, np.arange( i * f1o.NUM_basis, (i+1) * f1o.NUM_basis))
        np.testing.assert_array_equal(
            G2[i].full_vector, np.arange( i * f2.NUM_basis, (i+1) * f2.NUM_basis))

    if 1 in mesh.elements:
        np.testing.assert_array_equal(T1O[1].full_vector,
            np.array([3,4,5,10,11,12,13,14,15,16]))
    if 4 in mesh.elements:
        np.testing.assert_array_equal(T1O[4].full_vector,
            np.array([27,28,29,32,33,34,15,16,35,36]))

    f0 = FC('0-f-i', is_hybrid=False)
    f1i = FC('1-f-i', is_hybrid=False)
    f1o = FC('1-f-o', is_hybrid=False)
    f2 = FC('2-f-i', is_hybrid=False)

    G0 = f0.numbering.gathering
    G1I = f1i.numbering.gathering
    G1O = f1o.numbering.gathering
    G2 = f2.numbering.gathering

    if 1 in mesh.elements:
        np.testing.assert_array_equal(G0[1].full_vector,
            np.array([2, 12, 16, 5, 13, 17, 8, 14, 18, 11, 15, 19]))
    if 5 in mesh.elements:
        np.testing.assert_array_equal(G0[5].full_vector,
            np.array([19, 23, 27, 40, 43, 46, 41, 44, 47, 42, 45, 48]))

    if 1 in mesh.elements:
        np.testing.assert_array_equal(G1I[1].full_vector,
            np.array([17,21,18,22,19,23,20,24,10,25,28,13,26,29,16,27,30]))
    if 3 in mesh.elements:
        np.testing.assert_array_equal(G1I[3].full_vector,
            np.array([6,7,45,48,46,49,47,50,51,54,57,52,55,58,53,56,59]))

    if 2 in mesh.elements:
        np.testing.assert_array_equal(G1O[2].full_vector,
            np.array([20,31,34,21,32,35,22,33,36,37,41,38,42,39,43,40,44]))
    if 4 in mesh.elements:
        np.testing.assert_array_equal(G1O[4].full_vector,
            np.array([51,60,63,52,61,64,53,62,65,26,30,66,69,67,70,68,71]))

    for i in mesh.elements:
        np.testing.assert_array_equal(
            G2[i].full_vector, np.arange( i * f2.NUM_basis, (i+1) * f2.NUM_basis))

    return 1


def test_Form_NO3_mass_matrices():
    """"""
    if rAnk == mAster_rank:
        print("~~~ [test_Form_NO3_mass_matrices] ...... ", flush=True)

    mesh = MeshGenerator('crazy', c=0.12, bounds=((0,2), (0,3)))([2, 3], EDM='debug')
    space = SpaceInvoker('polynomials')([('Lobatto', 3), ('Lobatto', 2)])
    FC = FormCaller(mesh, space)

    f0 = FC('0-f-i', is_hybrid=False)
    f1i = FC('1-f-i', is_hybrid=False)
    f1o = FC('1-f-o', is_hybrid=False)
    f2 = FC('2-f-i', is_hybrid=False)

    M0 = f0.matrices.mass
    M1I = f1i.matrices.mass
    M1O = f1o.matrices.mass
    M2 = f2.matrices.mass

    if 0 in M1I:
        np.testing.assert_array_almost_equal(M1I[0].toarray()[0],
            np.array([ 0.60122847, -0.07187971,  0.01084695,  0.31006408, -0.03502005,
                       0.00256356, -0.12646503,  0.01871769, -0.00541012, -0.04149325,
                      -0.07118773,  0.01640982, -0.00061779,  0.01947493,  0.0279004 ,
                      -0.00217755, -0.00156191]))

        np.testing.assert_array_almost_equal(M1I[0].toarray()[5],
            np.array([  2.56355608e-03, -3.49248356e-02,  3.39292002e-01,  8.79317872e-02,
                       -3.43667975e-01,  3.31907618e+00,  1.90769236e-02, -4.86792311e-02,
                        4.53067783e-01, -5.71518131e-04,  1.25425544e-03,  1.22704151e-01,
                        1.25499757e-01,  6.05666173e-03, -3.39883326e-02,  2.98394338e-01,
                        2.64620426e-01]))
        np.testing.assert_array_almost_equal(M0[0].toarray()[4],
            np.array([  0.00554578,  0.00184104, -0.00249636,  0.00098171,  0.05000453,
                        0.02104327, -0.0189073 ,  0.006915  ,  0.00648875,  0.00315156,
                       -0.00204593,  0.00072488]))
        np.testing.assert_array_almost_equal(M2[0].toarray()[4],
            np.array([ 0.1989124 , -1.01161361,  0.17135838, -1.37435142,  7.3285936 , -1.53333516]))

    if 5 in M1I:
        np.testing.assert_array_almost_equal(M1I[5].toarray()[16],
            np.array([ -0.00246943,  0.00999065, -0.04423281,  0.01011988, -0.04317603,
                        0.15730403,  0.00479624, -0.02344232,  0.05624134, -0.00448805,
                        0.01015547, -0.00717787, -0.02391952,  0.03549893, -0.08172606,
                        0.04729042,  0.17342575]))
        np.testing.assert_array_almost_equal(M1I[5].toarray()[11],
            np.array([  0.00265893,  0.03943992,  0.20932935, -0.00753896,  0.09378854,
                        0.3155266 ,  0.00499139, -0.03020874, -0.0489097 , -0.07634322,
                        0.114424  ,  0.80176059,  0.06348914,  0.01112712, -0.02682325,
                       -0.13038107, -0.00717787]))
        np.testing.assert_array_almost_equal(M1O[5].toarray()[11],
            np.array([ 0.126342,  0.157021,  0.007539, -0.015498,  0.067897,  0.052683,
                       0.028214, -0.01012 ,  0.223465, -0.039418,  0.014051,  1.739309,
                      -0.30496 ,  0.129221,  0.233822, -0.036199,  0.01639 ]))

        np.testing.assert_array_almost_equal(M0[5].toarray()[11],
            np.array([ -0.00036701,  0.00063908, -0.00052651, -0.00175327,  0.0006056 ,
                       -0.00105295,  0.00170828,  0.00397803,  0.0013024 , -0.00237512,
                        0.00323713,  0.00803215]))
        np.testing.assert_array_almost_equal(M0[5].toarray()[7],
            np.array([  0.00086242, -0.00150339,  0.00039775,  0.00303505,  0.00578341,
                       -0.00948722,  0.00735124,  0.02618595,  0.0006056 , -0.00105295,
                        0.00170828,  0.00397803]))

    if 3 in M1I:
        np.testing.assert_array_almost_equal(M2[3].toarray()[0],
            np.array([11.5256843 , -1.47127256,  0.57304065, -1.58437845,  0.2166255 , -0.08919671]))
        np.testing.assert_array_almost_equal(M2[3].toarray()[1],
            np.array([-1.47127256,  5.64509819, -1.25453527,  0.2166255 , -0.9143819 , 0.2166255 ]))
        np.testing.assert_array_almost_equal(M2[3].toarray()[2],
            np.array([ 0.57304065, -1.25453527,  8.15018015, -0.08919671,  0.2166255 , -1.58437845]))
        np.testing.assert_array_almost_equal(M2[3].toarray()[3],
            np.array([-1.58437845,  0.2166255 , -0.08919671,  8.15018015, -1.25453527, 0.57304065]))
        np.testing.assert_array_almost_equal(M2[3].toarray()[4],
            np.array([ 0.2166255 , -0.9143819 ,  0.2166255 , -1.25453527,  5.64509819, -1.47127256]))
        np.testing.assert_array_almost_equal(M2[3].toarray()[5],
            np.array([-0.08919671,  0.2166255 , -1.58437845,  0.57304065, -1.47127256, 11.5256843 ]))

    return 1


def test_Form_NO4_Hodge():
    """"""
    if rAnk == mAster_rank:
        print("~~~ [test_Form_NO4_Hodge] ...... ", flush=True)

    def p(t, x, y): return np.cos(2*np.pi*x) * np.cos(2*np.pi*y) + t/2

    def P(t, x, y): return - np.pi * np.sin(2*np.pi*x) * np.cos(2*np.pi*y) + t/2
    def Q(t, x, y): return np.pi * np.cos(2*np.pi*x) * np.sin(2*np.pi*y) + t

    mesh = MeshGenerator('crazy', c=0.1, bounds=((0,1), (0,1)))([8, 9])
    space = SpaceInvoker('polynomials')([('Lobatto', 4), ('Lobatto', 3)])
    FC = FormCaller(mesh, space)

    scalar = FC('scalar', p)
    vector = FC('vector', (P, Q))
    t = random.uniform(-1, 1)

    f0 = FC('0-f-i', is_hybrid=False)
    f1i = FC('1-f-i', is_hybrid=False)
    f1o = FC('1-f-o', is_hybrid=False)
    f2 = FC('2-f-i', is_hybrid=False)

    # ------- f0 f2 -------------------------
    f0.TW.func.___DO_set_func_body_as___(scalar)
    f0.TW.___DO_push_all_to_instant___(t)
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
    f2.TW.func.___DO_set_func_body_as___(scalar)
    f2.TW.___DO_push_all_to_instant___(t)
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
        f0_cochain[i] = H @ f2.cochain.local[i]

    f0.cochain.local = f0_cochain
    f0_error = f0.error.L(upon=False)
    assert f0_error < 0.0015

    # ---------- f1 inner f1 outer ----------------
    f1i.TW.func.___DO_set_func_body_as___(vector)
    f1i.TW.___DO_push_all_to_instant___(t)
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
    f1o.TW.func.___DO_set_func_body_as___(vector)
    f1o.TW.___DO_push_all_to_instant___(t)
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
        f1i_cochain[i] = H @ f1o.cochain.local[i]

    f1i.cochain.local = f1i_cochain
    f1i_error = f1i.error.L(upon=False)
    assert f1i_error < 0.006
    return 1






if __name__ == '__main__':
    # mpiexec python _2dCSCG\TESTS\unittest_form.py
    # test_Form_NO4_Hodge()
    test_Form_NO1_coboundary()
    # test_Form_NO2_Naive_numbering()
    # test_Form_NO3_mass_matrices()