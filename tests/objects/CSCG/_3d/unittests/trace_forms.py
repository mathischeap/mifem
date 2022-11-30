# -*- coding: utf-8 -*-
import sys
if './' not in sys.path: sys.path.append('./')

from root.config.main import *
import random
from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller
from tests.objects.CSCG._3d.randObj.form_caller import random_FormCaller_of_total_load_around
# from scipy.sparse import csc_matrix
# from TOOLS.linear_algebra.data_structures import LocallyFullVector

def test_trace_NO__general_tests():
    """"""
    if RANK == MASTER_RANK:
        print("ttt [test_trace_NO__general_tests] ...... ", flush=True)

    # random form caller -----------------------------------------------------------------------------------
    if RANK == MASTER_RANK:
        load = random.randint(100, 499)
    else:
        load= None
    load = COMM.bcast(load, root=MASTER_RANK)
    FC = random_FormCaller_of_total_load_around(load, exclude_periodic=True)

    # --------- trace forms -------------------------------------------------------
    t0 = FC('0-t')
    t1 = FC('1-t')
    t2 = FC('2-t')

    # some tests ----------------------------------------------
    assert t0.whether.hybrid
    assert t1.whether.hybrid
    assert t2.whether.hybrid

    return 1



def test_trace_NO0_trace_0_form_Rd_and_Rc():
    """"""
    if RANK == MASTER_RANK:
        print("+0+ [test_trace_NO0_trace_0_form_Rd_and_Rc] ...... ", flush=True)

    if RANK == MASTER_RANK:
        load = random.randint(100, 499)
        t = random.random() * 10
    else:
        load= None
        t = None
    load, t = COMM.bcast([load, t], root=MASTER_RANK)

    #----------------- use crazy mesh ----------------------------------------
    FC = random_FormCaller_of_total_load_around(load, exclude_periodic=True)
    def pressure(t, x, y, z):
        return np.cos(1.5*np.pi*x) * \
               np.sin(2*np.pi*y) * \
               np.sin(np.pi*z-0.125)**2 * (1.25 - np.sin(t/2))
    P = FC('scalar', pressure)
    f0 = FC('0-f', is_hybrid=True)
    t0 = FC('0-t')
    f0.CF = P
    f0.CF.current_time = t
    t0.CF = P
    t0.CF.current_time = t

    f0.discretize() # default discretization, discrete a scalar to 0-form
    t0.discretize() # default discretization, discrete a scalar to 0-trace-form

    S = t0.space.selective_matrix._3dCSCG_0Trace[1]
    trace_map = t0.mesh.trace.elements.map
    for i in trace_map: # go through all local elements
        for j, side in enumerate('NSWEBF'): # go through all local trace elements around mesh element #i
            cochain_trace = t0.cochain.local_TEW[trace_map[i][j]]
            s = S[side]
            cochain_form = f0.cochain.local[i]
            cochain_trace_selective = s @ cochain_form
            np.testing.assert_array_almost_equal(cochain_trace, cochain_trace_selective)

    Sme = t0.matrices.selective # mesh-element-wise Selective matrix
    for i in t0.mesh.elements:
        Si = Sme[i]
        cochain_form = f0.cochain.local[i]
        cochain_trace_selective = Si @ cochain_form

        cochain_trace = t0.cochain.local[i]
        np.testing.assert_array_almost_equal(cochain_trace, cochain_trace_selective)

    return 1



def test_trace_NO1_trace_1_form_Rd_and_Rc():
    if RANK == MASTER_RANK:
        print("+1+ [test_trace_NO1_trace_1_form_Rd_and_Rc] ...... ", flush=True)

    def uuu(t, x, y, z):
        return np.cos(0.89*np.pi*x) + np.cos(np.pi*y) + np.cos(np.pi*z-0.578)**2 + np.sin(t)
    def vvv(t, x, y, z):
        return np.cos(2.21*np.pi*x) + np.cos(np.pi*y) + np.cos(np.pi*z-0.12)**2 * (1.5 + np.cos(t))
    def www(t, x, y, z):
        return np.cos(np.pi*x) + np.cos(np.pi*y) + np.cos(np.pi*z-0.15)**2 / (1.5 - np.sin(t))

    if RANK == MASTER_RANK:
        load = random.randint(100, 500)
        t = random.random() * 10
    else:
        load= None
        t = None
    load, t = COMM.bcast([load, t], root=MASTER_RANK)
    FC = random_FormCaller_of_total_load_around(load, exclude_periodic=True)
    velocity = FC('vector', (uuu, vvv, www))

    # first we test the Selective matrix --------------------------------------------------
    t1 = FC('1-t')
    t1.CF = velocity
    t1.CF.current_time = t
    t1.discretize() # Using the default T_para discretization

    f1 = FC('1-f', is_hybrid=True)
    f1.CF = velocity
    f1.CF.current_time = t
    f1.discretize()  # default discretization, discrete a vector to a standard 1-form

    S = t1.space.selective_matrix._3dCSCG_1Trace[1]

    trace_map = t1.mesh.trace.elements.map
    for i in trace_map: # go through all local elements
        for j, side in enumerate('NSWEBF'): # go through all local trace elements around mesh element #i
            cochain_trace = t1.cochain.local_TEW[trace_map[i][j]]
            s = S[side]
            cochain_form = f1.cochain.local[i]
            cochain_trace_selective = s @ cochain_form
            np.testing.assert_array_almost_equal(cochain_trace, cochain_trace_selective)

    Sme = t1.matrices.selective # mesh-element-wise Selective matrix
    for i in t1.mesh.elements: # go through all local mesh elements
        Si = Sme[i]
        cochain_form = f1.cochain.local[i]
        cochain_trace_selective = Si @ cochain_form

        cochain_trace = t1.cochain.local[i]
        np.testing.assert_array_almost_equal(cochain_trace, cochain_trace_selective)

    # now we compare that discretization from the vector is the same as discretization from its parallel component-------
    if RANK == MASTER_RANK:
        c = random.random() / 10
        if c < 0.05: c = 0
    else:
        c = None

    c = COMM.bcast(c, root=MASTER_RANK)
    mesh = MeshGenerator('crazy', c=c)([5, 4, 3], EDM='debug')
    space = SpaceInvoker('polynomials')([('Lobatto', 3), ('Lobatto', 4), ('Lobatto', 5)])
    FC = FormCaller(mesh, space)
    velocity = FC('vector', (uuu, vvv, www))
    t = 0

    t1 = FC('1-t')
    t1.CF = velocity
    t1.CF.current_time = t
    t1.discretize() # Using the default T_para discretization
    c1 = t1.cochain.local_TEW

    T1 = FC('1-t')
    para_V = velocity.components.T_para
    T1.CF = para_V
    T1.CF.current_time = t
    T1.discretize() # para_V is 'trace-element-wise', we use the trace-element-wise discretization
    C1 = T1.cochain.local_TEW

    for i in c1: assert np.max(np.abs(c1[i]-C1[i])) < 0.035

    return 1



def test_trace_NO2_trace_2_form_Rd_and_Rc():
    """"""
    if RANK == MASTER_RANK:
        print("+2+ [test_trace_NO2_trace_2_form_Rd_and_Rc] ...... ", flush=True)
    # test 2-Trace form, reduction from a standard vector...............

    if RANK == MASTER_RANK:
        load = random.randint(100, 499)
        t = random.random() * 10
    else:
        load= None
        t = None
    load, t = COMM.bcast([load, t], root=MASTER_RANK)

    #----------------- use crazy mesh ----------------------------------------
    FC = random_FormCaller_of_total_load_around(load, exclude_periodic=True)

    def u(t, x, y, z): return np.cos(np.pi*x) + np.sin(np.pi*y) * np.sin(np.pi*z-0.125)**2 + t/2
    def v(t, x, y, z): return np.sin(np.pi*x) + np.sin(np.pi*y) * np.sin(np.pi*z-0.125)**2 + t/2
    def w(t, x, y, z): return np.sin(np.pi*x) + np.cos(np.pi*y) * np.cos(np.pi*z-0.125)**2 + t
    vector = FC('vector', (u, v, w))
    f2 = FC('2-f', is_hybrid=True)
    t2 = FC('2-t')
    f2.CF = vector
    f2.CF.current_time = t
    t2.CF = vector
    t2.CF.current_time = t

    f2.discretize() # default discretization, discrete a vector to 2-form
    t2.discretize() # default discretization, discrete the outward norm component of a vector to 2-form

    S = t2.space.selective_matrix._3dCSCG_2Trace[1] # cannot use t2.matrices.selective

    trace_map = t2.mesh.trace.elements.map
    for i in trace_map: # go through all local mesh elements
        for j, side in enumerate('NSWEBF'): # go through all local trace elements around mesh element #i
            cochain_trace = t2.cochain.local_TEW[trace_map[i][j]]
            s = S[side]
            cochain_form = f2.cochain.local[i]
            cochain_trace_selective = s @ cochain_form
            np.testing.assert_array_almost_equal(cochain_trace, cochain_trace_selective)

    Sme = t2.matrices.selective # mesh-element-wise Selective matrix
    for i in t2.mesh.elements: # go through all local mesh elements
        Si = Sme[i]
        cochain_form = f2.cochain.local[i]
        cochain_trace_selective = Si @ cochain_form
        cochain_trace = t2.cochain.local[i]
        np.testing.assert_array_almost_equal(cochain_trace, cochain_trace_selective)

    return 1





def test_trace_NO3_trace_matrices():
    """"""
    if RANK == MASTER_RANK:
        load = random.randint(100, 299)
    else:
        load= None
    load = COMM.bcast(load, root=MASTER_RANK)
    FC = random_FormCaller_of_total_load_around(load, exclude_periodic=False)

    if RANK == MASTER_RANK:
        print(f"+++ [test_trace_NO3_trace_matrices] @load = {load} ...... ", flush=True)

    # Below we do a check for the assembled trace matrix for 2-trace-form---------------------------------
    t2 = FC('2-t')
    f2 = FC('2-f', is_hybrid=True)

    T2 = t2.matrices.trace
    T2.gathering_matrices = (t2, f2)
    T2a = T2.assembled
    M = T2a.do.gather_M_to_core(core=MASTER_RANK)
    SPL = f2.numbering.sharing_physical_locations
    if RANK == MASTER_RANK:
        COUNT = 0
        for Mi in M:
            NNZ = Mi.nnz
            if NNZ == 1:
                pass
            elif NNZ == 2:
                if tuple(Mi.indices) in SPL:
                    COUNT += 1
                else:
                    raise Exception()
            else:
                raise Exception()

        assert len(SPL) == COUNT


    # Below we do a check for the assembled trace matrix for 1-trace-form---------------------------------
    t1 = FC('1-t')
    f1 = FC('1-f', is_hybrid=True)

    T1 = t1.matrices.trace
    T1.gathering_matrices = (t1, f1)
    T1a = T1.assembled
    M = T1a.do.gather_M_to_core(core=MASTER_RANK)
    SPL = f1.numbering.sharing_physical_locations
    if RANK == MASTER_RANK:
        spl = list()
        for Mi in M:
            NNZ = Mi.nnz
            if NNZ == 1:
                pass
            elif NNZ == 2:
                i0, i1 = Mi.indices
                FOUND_IT = False
                for s in spl:
                    if i0 in s or i1 in s:
                        s.update((i0, i1))
                        FOUND_IT = True
                        break
                if FOUND_IT:
                    pass
                else:
                    spl.append({i0, i1})

            else:
                raise Exception()

        for s in spl:
            S = list(s)
            S.sort()
            assert tuple(S) in SPL, f"Something is wrong."

    # Below we do a check for the assembled trace matrix for 0-trace-form---------------------------------
    t0 = FC('0-t')
    f0 = FC('0-f', is_hybrid=True)

    T0 = t0.matrices.trace
    T0.gathering_matrices = (t0, f0)
    T0a = T0.assembled
    M = T0a.do.gather_M_to_core(core=MASTER_RANK)
    SPL = f0.numbering.sharing_physical_locations
    if RANK == MASTER_RANK:
        spl = list()
        for Mi in M:
            NNZ = Mi.nnz
            if NNZ == 1:
                pass
            elif NNZ == 2:
                i0, i1 = Mi.indices
                FOUND_IT = False
                for s in spl:
                    if i0 in s or i1 in s:
                        s.update((i0, i1))
                        FOUND_IT = True
                        break
                if FOUND_IT:
                    pass
                else:
                    spl.append({i0, i1})

            else:
                raise Exception()

        for s in spl:
            S = list(s)
            S.sort()
            assert tuple(S) in SPL, f"Something is wrong."

    return 1










if __name__ == '__main__':
    # mpiexec -n 6 python _3dCSCG\TESTS\unittest_trace.py
    # test_trace_NO__general_tests()
    test_trace_NO0_trace_0_form_Rd_and_Rc()
    test_trace_NO1_trace_1_form_Rd_and_Rc()
    test_trace_NO2_trace_2_form_Rd_and_Rc()
    # test_trace_NO3_trace_matrices()