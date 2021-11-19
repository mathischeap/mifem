
import sys
if './' not in sys.path: sys.path.append('./')

from root.config import *
import random
from _3dCSCG.main import MeshGenerator, SpaceInvoker, FormCaller#, ExactSolutionSelector


def test_trace_NO1_general_tests():
    """"""
    if rAnk == mAster_rank:
        print("ttt [test_trace_NO1_general_tests] ...... ", flush=True)

    mesh = MeshGenerator('crazy')([3, 3, 3], EDM=None, show_info=False)
    space = SpaceInvoker('polynomials')([('Lobatto', 2), ('Lobatto', 1), ('Lobatto', 2)], show_info=False)
    FC = FormCaller(mesh, space)

    t2 = FC('2-t')

    assert t2.IS_hybrid

    return 1


def test_trace_NO2_trace_form_Rd_and_Rc():
    """"""
    if rAnk == mAster_rank:
        print("ttt [test_trace_NO2_trace_form_Rd_and_Rc] ...... ",
              flush=True)
    # test 2-Trace form, reduction from a standard vector...............

    if rAnk == mAster_rank:
        el1 = random.randint(1,3)
        el2 = random.randint(1,3)
        el3 = random.randint(1,3)
        c = random.uniform(0.0, 0.3)
        if c < 0.1:c = 0
        t = random.uniform(0.0, 3)
    else:
        el1, el2, el3, c, t = [None for _ in range(5)]
    el1, el2, el3, c, t = cOmm.bcast([el1, el2, el3, c, t], root=mAster_rank)

    mesh = MeshGenerator('crazy', c=c)([el1, el2, el3])
    space = SpaceInvoker('polynomials')([('Lobatto', el3),
                                         ('Lobatto', el2),
                                         ('Lobatto', el1)])
    FC = FormCaller(mesh, space)
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

    S = space.selective_matrix._2Trace[1]

    trace_map = mesh.trace.elements.map
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
    test_trace_NO2_trace_form_Rd_and_Rc()