
import sys
if './' not in sys.path: sys.path.append('./')

from root.config import *
import random
from _3dCSCG.main import MeshGenerator, SpaceInvoker, FormCaller
from _3dCSCG.TESTS.random_objects import random_3D_FormCaller_of_total_load_around

def test_trace_NO__general_tests():
    """"""
    if rAnk == mAster_rank:
        print("ttt [test_trace_NO__general_tests] ...... ", flush=True)

    # random form caller -----------------------------------------------------------------------------------
    if rAnk == mAster_rank:
        load = random.randint(100, 499)
    else:
        load= None
    load = cOmm.bcast(load, root=mAster_rank)
    FC = random_3D_FormCaller_of_total_load_around(load, exclude_periodic=True)

    # --------- trace forms -------------------------------------------------------
    t0 = FC('0-t')
    t1 = FC('1-t')
    t2 = FC('2-t')

    # some tests ----------------------------------------------
    assert t0.IS_hybrid
    assert t1.IS_hybrid
    assert t2.IS_hybrid

    return 1



def test_trace_NO0_trace_0_form_Rd_and_Rc():
    """"""
    if rAnk == mAster_rank:
        print("+0+ [test_trace_NO0_trace_0_form_Rd_and_Rc] ...... ", flush=True)

    if rAnk == mAster_rank:
        load = random.randint(100, 499)
        t = random.random() * 10
    else:
        load= None
        t = None
    load, t = cOmm.bcast([load, t], root=mAster_rank)

    #----------------- use crazy mesh ----------------------------------------
    FC = random_3D_FormCaller_of_total_load_around(load, exclude_periodic=True)
    def pressure(t, x, y, z): return np.cos(1.5*np.pi*x) * np.sin(2*np.pi*y) * np.sin(np.pi*z-0.125)**2 * (1.25 - np.sin(t/2))
    P = FC('scalar', pressure)
    f0 = FC('0-f', is_hybrid=True)
    t0 = FC('0-t')
    f0.TW.func.DO.set_func_body_as(P)
    f0.TW.current_time = t
    f0.TW.DO.push_all_to_instant()
    t0.TW.func.DO.set_func_body_as(P)
    t0.TW.current_time = t
    t0.TW.DO.push_all_to_instant()

    f0.discretize() # default discretization, discrete a scalar to 0-form
    t0.discretize() # default discretization, discrete a scalar to 0-trace-form

    S = t0.space.selective_matrix._0Trace[1]
    trace_map = t0.mesh.trace.elements.map
    for i in trace_map: # go through all local elements
        for j, side in enumerate('NSWEBF'): # go through all local trace elements around mesh element #i
            cochain_trace = t0.cochain.local_TEW[trace_map[i][j]]
            s = S[side]
            cochain_form = f0.cochain.local[i]
            cochain_trace_selective = s @ cochain_form
            np.testing.assert_array_almost_equal(cochain_trace, cochain_trace_selective)

    return 1



def test_trace_NO1_trace_1_form_Rd_and_Rc():
    if rAnk == mAster_rank:
        print("+1+ [test_trace_NO1_trace_1_form_Rd_and_Rc] ...... ", flush=True)

    def uuu(t, x, y, z):
        return np.cos(0.89*np.pi*x) + np.cos(np.pi*y) + np.cos(np.pi*z-0.578)**2 + np.sin(t)
    def vvv(t, x, y, z):
        return np.cos(2.21*np.pi*x) + np.cos(np.pi*y) + np.cos(np.pi*z-0.12)**2 * (1.5 + np.cos(t))
    def www(t, x, y, z):
        return np.cos(np.pi*x) + np.cos(np.pi*y) + np.cos(np.pi*z-0.15)**2 / (1.5 - np.sin(t))

    if rAnk == mAster_rank:
        load = random.randint(100, 500)
        t = random.random() * 10
    else:
        load= None
        t = None
    load, t = cOmm.bcast([load, t], root=mAster_rank)
    FC = random_3D_FormCaller_of_total_load_around(load, exclude_periodic=True)
    velocity = FC('vector', (uuu, vvv, www))

    # first we test the Selective matrix --------------------------------------------------
    t1 = FC('1-t')
    t1.TW.func.DO.set_func_body_as(velocity)
    t1.TW.current_time = t
    t1.TW.DO.push_all_to_instant()
    t1.discretize() # Using the default T_para discretization

    f1 = FC('1-f', is_hybrid=True)
    f1.TW.func.DO.set_func_body_as(velocity)
    f1.TW.current_time = t
    f1.TW.DO.push_all_to_instant()
    f1.discretize()  # default discretization, discrete a vector to a standard 1-form

    S = t1.space.selective_matrix._1Trace[1]

    trace_map = t1.mesh.trace.elements.map
    for i in trace_map: # go through all local elements
        for j, side in enumerate('NSWEBF'): # go through all local trace elements around mesh element #i
            cochain_trace = t1.cochain.local_TEW[trace_map[i][j]]
            s = S[side]
            cochain_form = f1.cochain.local[i]
            cochain_trace_selective = s @ cochain_form
            np.testing.assert_array_almost_equal(cochain_trace, cochain_trace_selective)


    # now we compare that discretization from the vector is the same as discretization from its parallel component-------
    if rAnk == mAster_rank:
        c = random.random() / 10
        if c < 0.05: c = 0
    else:
        c = None

    c = cOmm.bcast(c, root=mAster_rank)
    mesh = MeshGenerator('crazy', c=c)([5, 4, 3], EDM='debug')
    space = SpaceInvoker('polynomials')([('Lobatto', 3), ('Lobatto', 4), ('Lobatto', 5)])
    FC = FormCaller(mesh, space)
    velocity = FC('vector', (uuu, vvv, www))
    t = 0

    t1 = FC('1-t')
    t1.TW.func.DO.set_func_body_as(velocity)
    t1.TW.current_time = t
    t1.TW.DO.push_all_to_instant()
    t1.discretize() # Using the default T_para discretization
    c1 = t1.cochain.local_TEW

    T1 = FC('1-t')
    para_V = velocity.components.T_para
    T1.TW.func.DO.set_func_body_as(para_V)
    T1.TW.current_time = t
    T1.TW.DO.push_all_to_instant()
    T1.discretize() # para_V is 'trace-element-wise', we use the trace-element-wise discretization
    C1 = T1.cochain.local_TEW

    for i in c1: assert np.max(np.abs(c1[i]-C1[i])) < 0.035

    return 1



def test_trace_NO2_trace_2_form_Rd_and_Rc():
    """"""
    if rAnk == mAster_rank:
        print("+2+ [test_trace_NO2_trace_2_form_Rd_and_Rc] ...... ", flush=True)
    # test 2-Trace form, reduction from a standard vector...............

    if rAnk == mAster_rank:
        load = random.randint(100, 499)
        t = random.random() * 10
    else:
        load= None
        t = None
    load, t = cOmm.bcast([load, t], root=mAster_rank)

    #----------------- use crazy mesh ----------------------------------------
    FC = random_3D_FormCaller_of_total_load_around(load, exclude_periodic=True)

    def u(t, x, y, z): return np.cos(np.pi*x) + np.sin(np.pi*y) * np.sin(np.pi*z-0.125)**2 + t/2
    def v(t, x, y, z): return np.sin(np.pi*x) + np.sin(np.pi*y) * np.sin(np.pi*z-0.125)**2 + t/2
    def w(t, x, y, z): return np.sin(np.pi*x) + np.cos(np.pi*y) * np.cos(np.pi*z-0.125)**2 + t
    vector = FC('vector', (u, v, w))
    f2 = FC('2-f', is_hybrid=True)
    t2 = FC('2-t')
    f2.TW.func.DO.set_func_body_as(vector)
    f2.TW.current_time = t
    f2.TW.DO.push_all_to_instant()
    t2.TW.func.DO.set_func_body_as(vector)
    t2.TW.current_time = t
    t2.TW.DO.push_all_to_instant()

    f2.discretize() # default discretization, discrete a vector to 2-form
    t2.discretize() # default discretization, discrete a the outward norm component of a vector to 2-form

    S = t2.space.selective_matrix._2Trace[1]

    trace_map = t2.mesh.trace.elements.map
    for i in trace_map: # go through all local elements
        for j, side in enumerate('NSWEBF'): # go through all local trace elements around mesh element #i
            cochain_trace = t2.cochain.local_TEW[trace_map[i][j]]
            s = S[side]
            cochain_form = f2.cochain.local[i]
            cochain_trace_selective = s @ cochain_form
            np.testing.assert_array_almost_equal(cochain_trace, cochain_trace_selective)

    return 1








if __name__ == '__main__':
    # mpiexec -n 12 python _3dCSCG\TESTS\unittest_trace.py
    # test_trace_NO__general_tests()
    # test_trace_NO0_trace_0_form_Rd_and_Rc()
    test_trace_NO1_trace_1_form_Rd_and_Rc()
    # test_trace_NO2_trace_2_form_Rd_and_Rc()