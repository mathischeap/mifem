import random
import sys
if './' not in sys.path: sys.path.append('./')

from root.config.main import *
from root.mifem.save import save, read
import os

from _3dCSCG.main import MeshGenerator, SpaceInvoker, FormCaller, ExactSolutionSelector
from _3dCSCG.tests.random_objects.form_caller import random_mesh_and_space_of_total_load_around


def test_ADF_NO1_general_tests_standard_forms():
    """"""
    if rAnk == mAster_rank:
        print("ADF [test_ADF_NO1_general_tests] ...... ", flush=True)

    mesh = MeshGenerator('crazy')([3, 3, 3], EDM=None, show_info=False)
    space = SpaceInvoker('polynomials')([2, 1, 2], show_info=False)
    FC = FormCaller(mesh, space)
    es = ExactSolutionSelector(mesh)('icpsNS:sincosRD')

    df0 = FC('0-adf', numbering_parameters={'scheme_name': 'Naive', })
    df1 = FC('1-adf', numbering_parameters={'scheme_name': 'Naive', })
    df2 = FC('2-adf', numbering_parameters={'scheme_name': 'Naive', })
    df3 = FC('3-adf', numbering_parameters={'scheme_name': 'Naive', })

    assert df0.IS.hybrid, "algebraic dual standard form must be hybrid."
    assert df1.IS.hybrid, "algebraic dual standard form must be hybrid."
    assert df2.IS.hybrid, "algebraic dual standard form must be hybrid."
    assert df3.IS.hybrid, "algebraic dual standard form must be hybrid."

    df0.prime.TW.func.do.set_func_body_as(es, 'pressure')
    df0.prime.TW.current_time = 0
    df0.prime.TW.do.push_all_to_instant()
    df0.prime.do.discretize()

    df0_error = df0.prime.error.L()
    save(df0, 'test_ADF_NO1_df0.mi')
    DF0 = read('test_ADF_NO1_df0.mi')
    DF0_error = DF0.prime.error.L()
    assert df0_error == DF0_error

    df1.prime.TW.func.do.set_func_body_as(es, 'vorticity')
    df1.prime.TW.current_time = 1
    df1.prime.TW.do.push_all_to_instant()
    df1.prime.do.discretize()
    df1_error = df1.prime.error.L()
    save(df1, 'test_ADF_NO1_df1.mi')
    DF1 = read('test_ADF_NO1_df1.mi')
    DF1_error = DF1.prime.error.L()
    assert df1_error == DF1_error

    df2.prime.TW.func.do.set_func_body_as(es, 'velocity')
    df2.prime.TW.current_time = 2
    df2.prime.TW.do.push_all_to_instant()
    df2.prime.do.discretize()
    df2_error = df2.prime.error.L()
    save(df2, 'test_ADF_NO1_df2.mi')
    DF2 = read('test_ADF_NO1_df2.mi')
    DF2_error = DF2.prime.error.L()
    assert df2_error == DF2_error

    df3.prime.TW.func.do.set_func_body_as(es, 'pressure')
    df3.prime.TW.current_time = 3
    df3.prime.TW.do.push_all_to_instant()
    df3.prime.do.discretize()
    df3_error = df3.prime.error.L()
    save(df3, 'test_ADF_NO1_df3.mi')
    DF3 = read('test_ADF_NO1_df3.mi')
    DF3_error = DF3.prime.error.L()
    assert df3_error == DF3_error


    if rAnk == mAster_rank:
        os.remove('test_ADF_NO1_df0.mi')
        os.remove('test_ADF_NO1_df1.mi')
        os.remove('test_ADF_NO1_df2.mi')
        os.remove('test_ADF_NO1_df3.mi')

    return 1




def test_ADF_NO2_general_tests_trace_forms():
    """"""
    if rAnk == mAster_rank:
        load = random.randint(100,1000)
        print(f"ADF [test_ADF_NO2_general_tests_trace_forms] @ load = {load}... ", flush=True)
    else:
        load = None
    load = cOmm.bcast(load, root=mAster_rank)
    mesh, space = random_mesh_and_space_of_total_load_around(load)
    FC = FormCaller(mesh, space)

    dt0 = FC('0-adt', numbering_parameters={'scheme_name': 'Naive', })
    dt1 = FC('1-adt', numbering_parameters={'scheme_name': 'Naive', })
    dt2 = FC('2-adt', numbering_parameters={'scheme_name': 'Naive', })

    assert dt0.IS.hybrid, "algebraic dual standard form must be hybrid."
    assert dt1.IS.hybrid, "algebraic dual standard form must be hybrid."
    assert dt2.IS.hybrid, "algebraic dual standard form must be hybrid."

    def u(t, x, y, z):
        return t + np.sin(2*np.pi*x) * np.cos(np.pi*y) * np.cos(2*np.pi*z)

    def v(t, x, y, z):
        return t + np.cos(2*np.pi*x) * np.sin(2*np.pi*y) * np.cos(np.pi*z)

    def w(t, x, y, z):
        return t + np.cos(np.pi*x) * np.cos(2*np.pi*y) * np.sin(2*np.pi*z)

    def p(t, x, y, z):
        return np.cos(2*np.pi*x) * np.cos(2*np.pi*y) * np.cos(2*np.pi*z) + t

    scalar = FC('scalar', p)
    vector = FC('vector', (u,v,w))

    # --------- save & read -----------------------------------------------
    dt0.prime.TW.func.do.set_func_body_as(scalar)
    dt0.prime.TW.current_time = 0
    dt0.prime.TW.do.push_all_to_instant()
    dt0.prime.do.discretize()
    dt1.prime.TW.func.do.set_func_body_as(vector)
    dt1.prime.TW.current_time = 0
    dt1.prime.TW.do.push_all_to_instant()
    dt1.prime.do.discretize()
    dt2.prime.TW.func.do.set_func_body_as(scalar)
    dt2.prime.TW.current_time = 0
    dt2.prime.TW.do.push_all_to_instant()
    dt2.prime.do.discretize()

    save([dt0, dt1, dt2], 'dual_objects.mi')

    DT0, DT1, DT2 = read('dual_objects.mi')

    assert DT0.prime.TW.func.body is None, f"when func is not a scalar or vector, it will not be saved."
    assert DT1.prime.TW.func.body is None, f"when func is not a scalar or vector, it will not be saved."
    assert DT2.prime.TW.func.body is None, f"when func is not a scalar or vector, it will not be saved."

    dt2c = dt2.cochain.local
    DT2c = DT2.cochain.local

    for i in dt2c:
        np.testing.assert_array_almost_equal(dt2c[i], DT2c[i])

    es = ExactSolutionSelector(mesh)('icpsNS:sincosRD')

    dt0.prime.TW.func.do.set_func_body_as(es, 'pressure')
    dt0.prime.TW.current_time = 0
    dt0.prime.TW.do.push_all_to_instant()
    dt0.prime.do.discretize()
    dt1.prime.TW.func.do.set_func_body_as(es, 'velocity')
    dt1.prime.TW.current_time = 2
    dt1.prime.TW.do.push_all_to_instant()
    dt1.prime.do.discretize()
    dt2.prime.TW.func.do.set_func_body_as(es, 'pressure')
    dt2.prime.TW.current_time = 3
    dt2.prime.TW.do.push_all_to_instant()
    dt2.prime.do.discretize()

    save([dt0, dt1, dt2], 'dual_objects.mi')
    Dt0, Dt1, Dt2 = read('dual_objects.mi', read_individuals=[1, 1, 1])

    if rAnk == mAster_rank:
        os.remove('dual_objects.mi')

    assert Dt0.prime.TW.func.body is not None
    assert Dt1.prime.TW.func.body is not None
    assert Dt2.prime.TW.func.body is not None

    dt2c = dt2.cochain.local
    Dt2c = Dt2.cochain.local

    for i in dt2c:
        np.testing.assert_array_almost_equal(dt2c[i], Dt2c[i])


    return 1





def test_ADF_NO3_coboundary():
    """"""
    if rAnk == mAster_rank:
        c = random.random() / 10
        if c < 0.05: c = 0
        print(f"ADF [test_ADF_NO3_coboundary] @ crazy mesh of c = %0.4f..."%c, flush=True)
    else:
        c = None
    c = cOmm.bcast(c, root=mAster_rank)
    mesh = MeshGenerator('crazy', c=c)([10, 9, 10], EDM=None, show_info=False)
    space = SpaceInvoker('polynomials')([3, 4, 3], show_info=False)
    FC = FormCaller(mesh, space)

    # ---- dual gradient ----------------------------------------------------------------------

    df3 = FC('3-adf', numbering_parameters={'scheme_name': 'Naive', })
    dt2 = FC('2-adt', numbering_parameters={'scheme_name': 'Naive', })
    def p(t, x, y, z): return np.sin(2*np.pi*x) * np.sin(2*np.pi*y) * np.sin(2*np.pi*z) + t
    def u(t, x, y, z): return 2*np.pi*np.cos(2*np.pi*x) * np.sin(2*np.pi*y) * np.sin(2*np.pi*z) + 0 * t
    def v(t, x, y, z): return 2*np.pi*np.sin(2*np.pi*x) * np.cos(2*np.pi*y) * np.sin(2*np.pi*z) + 0 * t
    def w(t, x, y, z): return 2*np.pi*np.sin(2*np.pi*x) * np.sin(2*np.pi*y) * np.cos(2*np.pi*z) + 0 * t

    scalar = FC('scalar', p)
    vector = FC('vector', (u, v, w))

    df3.prime.TW.func.do.set_func_body_as(scalar)
    df3.prime.TW.current_time = 0
    df3.prime.TW.do.push_all_to_instant()
    df3.prime.do.discretize()

    dt2.prime.TW.func.do.set_func_body_as(scalar)
    dt2.prime.TW.current_time = 0
    dt2.prime.TW.do.push_all_to_instant()
    dt2.prime.do.discretize()

    df2 = df3.coboundary(dt2)
    df2.prime.TW.func.do.set_func_body_as(vector)
    df2.prime.TW.current_time = 0
    df2.prime.TW.do.push_all_to_instant()

    assert df2.prime.error.L() < 0.01

    # ---- dual divergence ----------------------------------------------------------------------
    df1 = FC('1-adf')
    dt0 = FC('0-adt')

    def u(t, x, y, z): return np.cos(2*np.pi*x) * np.sin(2*np.pi*y) * np.sin(2*np.pi*z) + t
    def v(t, x, y, z): return np.sin(2*np.pi*x) * np.cos(2*np.pi*y) * np.sin(2*np.pi*z) + t
    def w(t, x, y, z): return np.sin(2*np.pi*x) * np.sin(2*np.pi*y) * np.cos(2*np.pi*z) + t
    def p(t, x, y, z): return - 6 * np.pi * np.sin(2*np.pi*x) * np.sin(2*np.pi*y) * np.sin(2*np.pi*z) + 0*t

    scalar = FC('scalar', p)
    vector = FC('vector', (u, v, w))

    df1.prime.TW.func.do.set_func_body_as(vector)
    df1.prime.TW.current_time = 1
    df1.prime.TW.do.push_all_to_instant()
    df1.prime.do.discretize()

    dt0.prime.TW.func.do.set_func_body_as(vector)
    dt0.prime.TW.current_time = 1
    dt0.prime.TW.do.push_all_to_instant()
    dt0.prime.do.discretize()

    df0 = df1.coboundary(dt0)
    df0.prime.TW.func.do.set_func_body_as(scalar)
    df0.prime.TW.current_time = 1
    df0.prime.TW.do.push_all_to_instant()

    assert df0.prime.error.L() < 0.028

    # ---- dual curl ----------------------------------------------------------------------
    mesh = MeshGenerator('crazy', c=0)([5, 6, 7], EDM=None, show_info=False)
    # currently, the coboundary for dual curl only works for orthogonal mesh.
    # The problem may be in the reduction, or mass matrix, or where-else.
    space = SpaceInvoker('polynomials')([4, 3, 2], show_info=False)
    FC = FormCaller(mesh, space)

    def F(t, x, y, z):
        return np.cos(np.pi*x) * np.sin(np.pi*y-0.125) * np.cos(np.pi*z) + t
    def G(t, x, y, z):
        return np.cos(np.pi*x) * np.cos(np.pi*y/2) * np.sin(np.pi*z) + t
    def H(t, x, y, z):
        return -2 * np.sin(np.pi*x) * np.cos(np.pi*y) * np.cos(np.pi*z) + t

    df2 = FC('2-adf')
    dt1 = FC('1-adt')

    V = FC('vector', (F, G, H))
    V_curl = V.numerical.curl
    V_perp = V.components.T_perp

    df2.prime.TW.func.do.set_func_body_as(V)
    df2.prime.TW.current_time = 0
    df2.prime.TW.do.push_all_to_instant()
    df2.prime.do.discretize()

    dt1.prime.TW.func.do.set_func_body_as(V_perp)
    dt1.prime.TW.current_time = 0
    dt1.prime.TW.do.push_all_to_instant()
    dt1.prime.do.discretize()

    df1 = df2.coboundary(dt1)
    df1.prime.TW.func.do.set_func_body_as(V_curl)
    df1.prime.TW.current_time = 0
    df1.prime.TW.do.push_all_to_instant()

    assert df2.prime.error.L() < 0.004, f"something is wrong."
    assert df1.prime.error.L() < 0.013, f"something is wrong."

    return 1









if __name__ == '__main__':
    # mpiexec -n 6 python _3dCSCG\tests\unittests\ADF.py
    # test_ADF_NO1_general_tests_standard_forms()
    test_ADF_NO2_general_tests_trace_forms()