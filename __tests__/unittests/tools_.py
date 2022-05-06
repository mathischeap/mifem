# -*- coding: utf-8 -*-
"""
@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
import sys
if './' not in sys.path: sys.path.append('./')

from root.config.main import *
import os
from time import sleep
from scipy import sparse as spspa
from tools.iterators.simple import SimpleIterator
import random
from tools.linear_algebra.data_structures.global_matrix.main import GlobalVector, DistributedVector, GlobalMatrix
from tools.linear_algebra.deprecated.data_structure_old0 import GlobalMatrix as GlobalMatrixOld
import tools.linear_algebra.deprecated.gmres as gmres
import tools.linear_algebra.solvers.serial.deprecated as serial_spspalinalg
from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller
from objects.CSCG._2d.master import MeshGenerator as MeshGenerator2D
from objects.CSCG._2d.master import SpaceInvoker as SpaceInvoker2D
from objects.CSCG._2d.master import FormCaller as FormCaller2D
from tools.linear_algebra.gathering.regular.chain_matrix.main import Chain_Gathering_Matrix
from tools.linear_algebra.elementwise_cache.operators.concatenate.main import bmat, concatenate
import tools.linear_algebra.elementwise_cache.operators.concatenate.main as mif
from tools.linear_algebra.linear_system.main import LinearSystem
from tools.linear_algebra.elementwise_cache.objects.sparse_matrix.main import EWC_ColumnVector
from tools.run.reader import ParallelMatrix3dInputRunner, RunnerDataReader
from objects.CSCG.base.forms.base.BC.partial_cochain.partial_dofs.main import PartialDofs
from objects.CSCG._3d.__tests__.Random.form_caller import random_mesh_and_space_of_total_load_around

def ___TEST_SOLVER___(tk, tk1):
    """
    Parameters
    ----------
    tk :
    tk1 ：

    Returns
    -------
    exit_code : The standard exit code.
    shut_down : If it is ``True``, the outer iterator will shutdown immediately.
    message : The solver message.
    output1 :
    output2 :
    """
    sleep(0.01)
    assert tk1 > tk
    if rAnk == mAster_rank:
        ot1 = (tk+tk1)/2
        ot2 = (tk+0.1*tk1)/3
    else:
        ot1 = None
        ot2 = None
    ot1, ot2 = cOmm.bcast([ot1, ot2], root=mAster_rank)
    return 1, 0, 'random message', ot1, ot2

def test_TOOLS_NO1_iterator():
    if rAnk == mAster_rank:
        print("=== [test_TOOLS_NO1_iterator] ...... ", flush=True)
    SI = SimpleIterator(t0=0, dt=0.1, max_steps=10,
                        auto_save_frequency=5,
                        monitor_factor=0,
                        RDF_filename='RDF_filename',
                        save_to_mitr=True,
                        name=None)
    SI(___TEST_SOLVER___, [0,0])
    SI.run()

    S2 = SimpleIterator.read("RDF_filename.mitr")

    if rAnk == mAster_rank:
        assert S2._solver_dir_ == SI._solver_dir_
        assert S2._solver_source_code_ == SI._solver_source_code_
        np.testing.assert_array_equal(SI.RDF.to_numpy(), S2.RDF.to_numpy())

        os.remove('RDF_filename.csv')
        os.remove('RDF_filename.mitr')
    else:
        assert S2 is None, "Iterator only read to master core. It is only for retrieve info, not for resume simulation."
    return 1



def test_TOOLS_NO2_0_linear_algebra_gmres0_solver():
    if rAnk == mAster_rank:
        print("=== [test_TOOLS_NO2_0_linear_algebra_gmres0_solver] ...... ", flush=True)
    A = np.array([(1, 4, 7),
                  (2, 9, 7),
                  (5, 8, 3)])
    b = np.array([(1,),
                  (8,),
                  (2,)])
    Ar = np.random.rand(3,3)
    br = np.random.rand(3,1)
    i = random.randint(0, 2)
    j = random.randint(0, 2)
    k = random.randint(0, 2)
    Ar[:,j] = 0
    Ar[i,:] = 0
    br[k] = 0
    if sIze > 3:
        if rAnk == sIze -1:
            Ar = np.zeros((3,3)) # A even can be empty in some cores.
    AA = cOmm.gather(Ar, root=0)
    bb = cOmm.gather(br, root=0)
    if rAnk == 0:
        a0 = np.zeros((3,3))
        b0 = np.zeros((3,1))
        for i in range(sIze):
            if i != 0:
                a0 += AA[i]
                b0 += bb[i]
        Ar = A - a0
        br = b - b0
    A = GlobalMatrixOld(spspa.csc_matrix(Ar))
    b = GlobalVector(spspa.csc_matrix(br))
    x0 = DistributedVector(spspa.csc_matrix((3,1)))
    x0, info = gmres.gmres0(A, b, x0, restart=3, tol=1e-6, maxiter=1)[0:2]
    if 0 in x0.indices: np.testing.assert_almost_equal(x0[0,0], -2.1810344827586)
    if 1 in x0.indices: np.testing.assert_almost_equal(x0[1,0],  1.8362068965517)
    if 2 in x0.indices: np.testing.assert_almost_equal(x0[2,0], -0.5948275862068)
    assert info == 0
    A.DO.___PRIVATE_being_regularly_distributed___('row') # this will gathering all data to core 0.
    if rAnk == 0:
        np.testing.assert_array_almost_equal(A.M.toarray(),
                                      np.array([(1, 4, 7),
                                                (2, 9, 7),
                                                (5, 8, 3)]))
    A = GlobalMatrixOld(spspa.csc_matrix(Ar))
    A.DO.___PRIVATE_being_regularly_distributed___('column') # this will gathering all data to core 0.
    if rAnk == 0:
        np.testing.assert_array_almost_equal(A.M.toarray(),
                                      np.array([(1, 4, 7),
                                                (2, 9, 7),
                                                (5, 8, 3)]))
    return 1

def test_TOOLS_NO2_1_linear_algebra_gmres1_solver():
    if rAnk == mAster_rank:
        print("=== [test_TOOLS_NO2_1_linear_algebra_gmres1_solver] ...... ", flush=True)
    A = np.array([(1, 4, 7),
                  (2, 9, 7),
                  (5, 8, 3)])
    b = np.array([(1,),
                  (8,),
                  (2,)])
    Ar = np.random.rand(3,3)
    br = np.random.rand(3,1)
    i = random.randint(0, 2)
    j = random.randint(0, 2)
    k = random.randint(0, 2)
    Ar[:,j] = 0
    Ar[i,:] = 0
    br[k] = 0
    if sIze > 3:
        if rAnk == sIze -1:
            Ar = np.zeros((3,3)) # A even can be empty in some cores.
    AA = cOmm.gather(Ar, root=0)
    bb = cOmm.gather(br, root=0)
    if rAnk == 0:
        a0 = np.zeros((3,3))
        b0 = np.zeros((3,1))
        for i in range(sIze):
            if i != 0:
                a0 += AA[i]
                b0 += bb[i]
        Ar = A - a0
        br = b - b0
    A = GlobalMatrixOld(spspa.csc_matrix(Ar))
    b = GlobalVector(spspa.csc_matrix(br))
    x0 = DistributedVector(spspa.csc_matrix((3,1)))

    x0, info = gmres.gmres1(A, b, x0, restart=3, tol=1e-6, maxiter=1)[0:2]

    x0 = x0.V
    if 0 in x0.indices: np.testing.assert_almost_equal(x0[0,0], -2.1810344827586)
    if 1 in x0.indices: np.testing.assert_almost_equal(x0[1,0],  1.8362068965517)
    if 2 in x0.indices: np.testing.assert_almost_equal(x0[2,0], -0.5948275862068)
    assert info == 0

    return 1

def test_TOOLS_NO2_2_linear_algebra_gmres2_solver():
    if rAnk == mAster_rank:
        print("=== [test_TOOLS_NO2_2_linear_algebra_gmres2_solver] ...... ", flush=True)
    A = np.array([(1, 4, 7),
                  (2, 9, 7),
                  (5, 8, 3)])
    b = np.array([(1,),
                  (8,),
                  (2,)])
    _A = A
    Ar = np.random.rand(3,3)
    br = np.random.rand(3,1)
    i = random.randint(0, 2)
    j = random.randint(0, 2)
    Ar[:,j] = 0
    Ar[i,:] = 0
    br[i,0] = 0
    if sIze > 2:
        if rAnk == sIze - 1:
            if i == 1:
                # sometimes, even one A can be empty in the last core.
                Ar = np.zeros((3,3))
                br[:,0] = 0
    AA = cOmm.gather(Ar, root=0)
    bb = cOmm.gather(br, root=0)
    if rAnk == 0:
        a0 = np.zeros((3,3))
        b0 = np.zeros((3,1))
        for i in range(sIze):
            if i != 0:
                a0 += AA[i]
                b0 += bb[i]
        Ar = A - a0
        br = b - b0

    if rAnk == 0:
        k = random.randint(0, 2)
        # k = 0
        Ar[k,:] = A[k,:]
        br[k,0] = b[k,0]
    else:
        k = None
    k = cOmm.bcast(k, root=0)
    if rAnk != 0:
        Ar[k,:] = 0
        br[k,:] = 0

    AAA = cOmm.gather(Ar, root=mAster_rank)
    if rAnk == mAster_rank: # check that the GMA indeed is a representation of A.
        AAA = np.sum(AAA, axis=0)
        np.testing.assert_array_almost_equal(AAA, A)

    A = GlobalMatrixOld(spspa.csc_matrix(Ar))
    b = GlobalVector(spspa.csc_matrix(br))
    x0 = DistributedVector(spspa.csc_matrix((3,1)))

    AAA = A.___PRIVATE_gather_M_to_core___()
    if rAnk == mAster_rank: # check that the GMA indeed is a representation of A.
        np.testing.assert_array_almost_equal(AAA.toarray(), _A)

    x0, info = gmres.gmres2(A, b, x0, restart=3, tol=1e-6, maxiter=1)[0:2]

    if 0 in x0.indices: np.testing.assert_almost_equal(x0[0,0], -2.1810344827586)
    if 1 in x0.indices: np.testing.assert_almost_equal(x0[1,0],  1.8362068965517)
    if 2 in x0.indices: np.testing.assert_almost_equal(x0[2,0], -0.5948275862068)
    assert info == 0

    return 1

def test_TOOLS_NO2_3_linear_algebra_serial_scipy_sparse_solver():
    if rAnk == mAster_rank:
        print("=== [test_TOOLS_NO2_3_linear_algebra_serial_scipy_sparse_solver] ...... ", flush=True)

    A = np.array([(1, 4, 7),
                  (2, 9, 7),
                  (5, 8, 3)])
    b = np.array([(1,),
                  (8,),
                  (2,)])
    Ar = np.random.rand(3,3)
    br = np.random.rand(3,1)
    i = random.randint(0, 2)
    j = random.randint(0, 2)
    k = random.randint(0, 2)
    Ar[:,j] = 0
    Ar[i,:] = 0
    br[k] = 0
    if sIze > 3:
        if rAnk == sIze -1:
            if i in (1,2):
                Ar = np.zeros((3,3)) # A even can be empty in some cores.
    AA = cOmm.gather(Ar, root=0)
    bb = cOmm.gather(br, root=0)
    if rAnk == 0:
        a0 = np.zeros((3,3))
        b0 = np.zeros((3,1))
        for i in range(sIze):
            if i != 0:
                a0 += AA[i]
                b0 += bb[i]
        Ar = A - a0
        br = b - b0
    A = GlobalMatrixOld(spspa.csc_matrix(Ar))
    b = GlobalVector(spspa.csc_matrix(br))
    X0 = DistributedVector(spspa.csc_matrix((3,1)))

    A = A.___PRIVATE_gather_M_to_core___(clean_local=True)
    A = GlobalMatrixOld(A)

    b = b.___PRIVATE_gather_V_to_core___(clean_local=True)
    b = GlobalVector(b)

    x0, info = getattr(serial_spspalinalg, 'gmres')(A, b, X0, restart=3, tol=1e-6, maxiter=1)[0:2]
    assert x0.IS.master_dominating, "Results must be a master dominating DistributedVector."
    x0 = x0.V
    if rAnk == mAster_rank:
        np.testing.assert_almost_equal(x0[0,0], -2.1810344827586)
        np.testing.assert_almost_equal(x0[1,0],  1.8362068965517)
        np.testing.assert_almost_equal(x0[2,0], -0.5948275862068)
    else:
        assert len(x0.indices) == 0, "solutions are in core master only."
    assert info == 0

    x0, info = getattr(serial_spspalinalg, 'gcrotmk')(A, b, X0, restart=3, tol=1e-6, maxiter=5)[0:2]
    assert x0.IS.master_dominating, "Results must be a master dominating DistributedVector."
    x0 = x0.V
    if rAnk == mAster_rank:
        np.testing.assert_almost_equal(x0[0,0], -2.1810344827586)
        np.testing.assert_almost_equal(x0[1,0],  1.8362068965517)
        np.testing.assert_almost_equal(x0[2,0], -0.5948275862068)
    else:
        assert len(x0.indices) == 0, "solutions are in core master only."
    assert info == 0

    x0, info = getattr(serial_spspalinalg, 'bicgstab')(A, b, X0, restart=3, tol=1e-6, maxiter=5)[0:2]
    assert x0.IS.master_dominating, "Results must be a master dominating DistributedVector."
    x0 = x0.V
    if rAnk == mAster_rank:
        np.testing.assert_almost_equal(x0[0,0], -2.1810344827586)
        np.testing.assert_almost_equal(x0[1,0],  1.8362068965517)
        np.testing.assert_almost_equal(x0[2,0], -0.5948275862068)
    else:
        assert len(x0.indices) == 0, "solutions are in core master only."
    assert info == 0

    return 1

def test_TOOLS_NO2_4_linear_algebra_serial_spsolve_solver():
    if rAnk == mAster_rank:
        print(">>> [test_TOOLS_NO2_4_linear_algebra_serial_spsolve_solver] ...... ", flush=True)

    A = np.array([(1, 4, 7),
                  (2, 9, 7),
                  (5, 8, 3)])
    b = np.array([(1,),
                  (8,),
                  (2,)])
    Ar = np.random.rand(3,3)
    br = np.random.rand(3,1)
    i = random.randint(0, 2)
    j = random.randint(0, 2)
    k = random.randint(0, 2)
    Ar[:,j] = 0
    Ar[i,:] = 0
    br[k] = 0
    if sIze > 3:
        if rAnk == sIze -1:
            if i in (1,2):
                Ar = np.zeros((3,3)) # A even can be empty in some cores.
    AA = cOmm.gather(Ar, root=0)
    bb = cOmm.gather(br, root=0)
    if rAnk == 0:
        a0 = np.zeros((3,3))
        b0 = np.zeros((3,1))
        for i in range(sIze):
            if i != 0:
                a0 += AA[i]
                b0 += bb[i]
        Ar = A - a0
        br = b - b0
    A = GlobalMatrixOld(spspa.csc_matrix(Ar))
    b = GlobalVector(spspa.csc_matrix(br))

    A = A.___PRIVATE_gather_M_to_core___(clean_local=True)
    A = GlobalMatrixOld(A)

    b = b.___PRIVATE_gather_V_to_core___(clean_local=True)
    b = GlobalVector(b)

    x0, info = serial_spspalinalg.spsolve(A, b)[0:2]

    assert x0.IS.master_dominating, "Results must be a master dominating DistributedVector."
    x0 = x0.V
    if rAnk == mAster_rank:
        np.testing.assert_almost_equal(x0[0,0], -2.1810344827586)
        np.testing.assert_almost_equal(x0[1,0],  1.8362068965517)
        np.testing.assert_almost_equal(x0[2,0], -0.5948275862068)
    else:
        assert len(x0.indices) == 0, "solutions are in core master only."
    assert info == 0

    return 1


def test_TOOLS_NO3_GlobalMatrix_GlobalVector_operators_test():
    if rAnk == mAster_rank:
        print("=== [test_TOOLS_NO3_GlobalMatrix_GlobalVector_operators_test] ....", flush=True)
    # we test C = A @ B, e = C @ d
    if rAnk == mAster_rank:
        i = random.randint(25, 49)
        j = random.randint(50, 74)
        k = random.randint(75, 99)
    else:
        i, j, k = None, None, None
    i = cOmm.bcast(i, root=mAster_rank)
    j = cOmm.bcast(j, root=mAster_rank)
    k = cOmm.bcast(k, root=mAster_rank)
    sa = random.randint(10, 100)/5000
    sb = random.randint(10, 100)/5000
    sd = random.randint(80, 100)/120
    A0 = spspa.random(i, j, sa, format='csr')
    A1 = A0.copy()
    A2 = A0.copy()
    A3 = A0.copy()
    B0 = spspa.random(j, k, sb, format='csc')
    B1 = B0.copy()
    B2 = B0.copy()
    B3 = B0.copy()
    d = spspa.random(k, 1, sd, format='csc')
    A0 = cOmm.gather(A0, root=mAster_rank)
    B0 = cOmm.gather(B0, root=mAster_rank)
    d0 = cOmm.gather(d, root=mAster_rank)
    if rAnk == mAster_rank:
        A0 = np.sum(A0)
        B0 = np.sum(B0)
        d0 = np.sum(d0)
        C0 = A0 @ B0
        e0 = (C0 @ d0).toarray()[:,0]
    del A0, B0
    A1 = GlobalMatrixOld(A1)
    A2 = GlobalMatrixOld(A2)
    A3 = GlobalMatrixOld(A3)
    B1 = GlobalMatrixOld(B1)
    B2 = GlobalMatrixOld(B2)
    B3 = GlobalMatrixOld(B3)

    assert A1.M.nnz == A2.M.nnz == A3.M.nnz
    assert B1.M.nnz == B2.M.nnz == B3.M.nnz
    assert A1.IS_regularly_distributed is False
    assert B1.IS_regularly_distributed is False
    C1 = A1 @ B1
    B2.DO.___PRIVATE_being_regularly_distributed___('column')
    assert A2.IS_regularly_distributed is False
    assert B2.IS_regularly_distributed == 'column'
    C2 = A2 @ B2

    A3.DO.___PRIVATE_being_regularly_distributed___('row')
    assert A3.IS_regularly_distributed == 'row'
    assert B3.IS_regularly_distributed is False
    C3 = A3 @ B3

    CM123 = 3 * C1 + C2 * 0.7 - C3 *2.1
    CM123 = CM123.___PRIVATE_gather_M_to_core___()

    M1 = C1.___PRIVATE_gather_M_to_core___()
    M2 = C2.___PRIVATE_gather_M_to_core___()
    M3 = C3.___PRIVATE_gather_M_to_core___()
    M1T = C1.T.___PRIVATE_gather_M_to_core___()

    if rAnk == mAster_rank:
        # noinspection PyUnboundLocalVariable
        np.testing.assert_almost_equal(np.sum(np.abs((3+0.7-2.1)*C0 - CM123)), 0)
        np.testing.assert_almost_equal(np.sum(np.abs(C0 - M1)), 0)
        np.testing.assert_almost_equal(np.sum(np.abs(C0 - M2)), 0)
        np.testing.assert_almost_equal(np.sum(np.abs(C0 - M3)), 0)
        np.testing.assert_almost_equal(np.sum(np.abs(C0.T - M1T)), 0)

    d = GlobalVector(d)

    e1 = C1 @ d
    e2 = C2 @ d
    e3 = C3 @ d
    E1 = e1.___PRIVATE_gather_V_to_core___()
    E2 = e2.___PRIVATE_gather_V_to_core___()
    E3 = e3.___PRIVATE_gather_V_to_core___()
    if rAnk == mAster_rank:
        # noinspection PyUnboundLocalVariable
        np.testing.assert_almost_equal(np.sum(np.abs(E1-e0)), 0)
        np.testing.assert_almost_equal(np.sum(np.abs(E2-e0)), 0)
        np.testing.assert_almost_equal(np.sum(np.abs(E3-e0)), 0)
    e12 = e1-e2
    e13 = e1-e3
    E12 = e12.___PRIVATE_gather_V_to_core___()
    E13 = e13.___PRIVATE_gather_V_to_core___()
    if rAnk == mAster_rank:
        np.testing.assert_almost_equal(np.sum(np.abs(E12)), 0)
        np.testing.assert_almost_equal(np.sum(np.abs(E13)), 0)
    eP123 = e1 + e2 + e3
    EP123 = eP123.___PRIVATE_gather_V_to_core___()
    if rAnk == mAster_rank:
        np.testing.assert_almost_equal(np.sum(np.abs(EP123-3*e0)), 0)

    eM123 = 0.5*e1 + e2*2 - 0.125*e3
    eM123 = eM123.___PRIVATE_gather_V_to_core___()
    if rAnk == mAster_rank:
        np.testing.assert_almost_equal(np.sum(np.abs(eM123-(0.5+2-0.125)*e0)), 0)

    return 1


def test_TOOLS_NO4_GlobalMatrix_GlobalVector_operators_test():
    """"""
    if rAnk == mAster_rank:
        print("=== [test_TOOLS_NO4_GlobalMatrix_GlobalVector_operators_test] ....", flush=True)
    # we test C = A @ B
    if rAnk == mAster_rank:
        i, j, k = random.randint(2, 16), random.randint(11, 20), random.randint(21, 30)
    else:
        i, j, k = None, None, None
    i, j, k = cOmm.bcast(i, root=mAster_rank), cOmm.bcast(j, root=mAster_rank), cOmm.bcast(k, root=mAster_rank)
    if rAnk == mAster_rank:
        rand_list = random.sample(range(0, i), i)
        C = [i // sIze + (1 if x < i % sIze else 0) for x in range(sIze)]
        not_empty = list()
        for c in C:
            eci = list()
            for d in range(c):
                eci.append(rand_list.pop())
            not_empty.append(eci)
    else:
        not_empty = None
    not_empty = cOmm.scatter(not_empty, root=mAster_rank)

    # column-major A, random B
    _A_ = spspa.random(j, i, random.randint(80, 100)/120, format='csc') # base A
    A = spspa.lil_matrix((j,i))
    A[:, not_empty] = _A_[:, not_empty]
    A = A.tocsc()
    B = spspa.random(i, k, random.randint(60, 80)/120)
    A0 = cOmm.gather(A, root=mAster_rank)
    B0 = cOmm.gather(B, root=mAster_rank)
    if rAnk == mAster_rank:
        A0 = np.sum(A0)
        B0 = np.sum(B0)
        C0 = A0 @ B0
    A = GlobalMatrixOld(A)
    A.IS_regularly_distributed = 'column'
    assert A.DO.___PRIVATE_check_if_Iam_column_major___() == (True, 0)
    B = GlobalMatrixOld(B)
    assert (A.IS_regularly_distributed, B.IS_regularly_distributed) == ('column', False)
    C = A @ B
    C = C.___PRIVATE_gather_M_to_core___()
    if rAnk == mAster_rank:
        # noinspection PyUnboundLocalVariable
        np.testing.assert_almost_equal(np.sum(np.abs(C - C0)), 0)

    # random B, row-major A: B @ A
    _A_ = spspa.random(i, j, random.randint(80, 100)/120, format='csr') # base A
    A = spspa.lil_matrix((i,j))
    A[not_empty, :] = _A_[not_empty, :]
    A = A.tocsr()
    B = spspa.random(k, i, random.randint(60, 80)/120)
    A0 = cOmm.gather(A, root=mAster_rank)
    B0 = cOmm.gather(B, root=mAster_rank)
    if rAnk == mAster_rank:
        A0 = np.sum(A0)
        B0 = np.sum(B0)
        C0 = B0 @ A0
    A = GlobalMatrixOld(A)
    A.IS_regularly_distributed = 'row'
    assert A.DO.___PRIVATE_check_if_Iam_row_major___() == (True, 0)
    B = GlobalMatrixOld(B)
    assert (B.IS_regularly_distributed, A.IS_regularly_distributed) == (False, 'row')
    C = B @ A
    C = C.___PRIVATE_gather_M_to_core___()
    if rAnk == mAster_rank:
        np.testing.assert_almost_equal(np.sum(np.abs(C - C0)), 0)

    # test do.claim_distribution_pattern
    if sIze == 1: # when single core, A is full, so it is to be considered as, always, row-major
        pass
    else:
        _A_ = spspa.random(i, i, 1, format='csc') # base A
        A = spspa.lil_matrix((i,i))
        A[:, not_empty] = _A_[:, not_empty]
        A = A.tocsc()
        A = GlobalMatrixOld(A)
        assert A.DO.claim_distribution_pattern() == ('column', 0)
        assert A.IS_regularly_distributed == 'column'

    _A_ = spspa.random(i, i, 1, format='csr')  # base A
    A = spspa.lil_matrix((i, i))
    A[not_empty,:] = _A_[not_empty,:]
    A = A.tocsr()
    A = GlobalMatrixOld(A)
    assert A.DO.claim_distribution_pattern() == ('row', 0)
    assert A.IS_regularly_distributed == 'row'

    # random: A @ B, (A.IS_regularly_distributed, B.IS_regularly_distributed) = (False, False)
    for _ in range(4): # test four times to test some different cases.
        if rAnk == mAster_rank:
            i = random.randint(2, 16)
        else:
            i = None
        i = cOmm.bcast(i, root=mAster_rank)
        A = spspa.random(i, i, random.randint(10, 50)/100, format='csr')
        B = spspa.random(i, i, random.randint(10, 50)/100, format='csc')
        A0 = cOmm.gather(A, root=mAster_rank)
        B0 = cOmm.gather(B, root=mAster_rank)
        if rAnk == mAster_rank:
            A0 = np.sum(A0)
            B0 = np.sum(B0)
            C0 = A0 @ B0
        A = GlobalMatrixOld(A)
        B = GlobalMatrixOld(B)
        assert (B.IS_regularly_distributed, A.IS_regularly_distributed) == (False, False)
        C = A @ B
        C = C.___PRIVATE_gather_M_to_core___()
        if rAnk == mAster_rank:
            np.testing.assert_almost_equal(np.sum(np.abs(C - C0)), 0)

    return 1



def test_TOOLS_NO5_DistributedVector_operators_test():
    """Unittest for iterators."""
    if rAnk == mAster_rank:
        print("=== [test_TOOLS_NO5_DistributedVector_operators_test] ....", flush=True)
    # y = A @ x
    if rAnk == mAster_rank:
        i = random.randint(25, 49)
        j = random.randint(50, 74)
        x0 = spspa.random(j, 1, 0.8, format='lil')
    else:
        x0 = None
        i, j = None, None
    i = cOmm.bcast(i, root=mAster_rank)
    j = cOmm.bcast(j, root=mAster_rank)
    sa = random.randint(10, 100)/500
    A = spspa.random(i, j, sa, format='csc')
    A0 = cOmm.gather(A, root=mAster_rank)
    if rAnk == mAster_rank:
        A0 = np.sum(A0)
        y0 = (A0 @ x0).toarray()[:,0]

    non_empty_list = [j // sIze + (1 if x < j % sIze else 0) for x in range(sIze)]
    x0 = cOmm.bcast(x0, root=mAster_rank)
    x = spspa.lil_matrix((j, 1))

    non_empty_range = range(int(np.sum(non_empty_list[:rAnk])), int(np.sum(non_empty_list[:rAnk+1])))
    non_empty_range = [i for i in non_empty_range]

    data = x0[non_empty_range, 0].toarray().ravel()
    x[non_empty_range, 0] = data

    A = GlobalMatrixOld(A)
    x = DistributedVector(x.tocsc())
    y1 = A @ x
    y1 = y1.___PRIVATE_gather_V_to_core___()
    if rAnk == mAster_rank:
        # noinspection PyUnboundLocalVariable
        np.testing.assert_array_almost_equal(y1-y0, 0)

    return 1

def test_TOOLS_NO6_send_GM_in_parts_test():
    if rAnk == mAster_rank:
        print("=== [test_TOOLS_NO6_send_GM_in_parts_test] ....", flush=True)

    for T in range(50):
        if rAnk == mAster_rank:
            i = random.randint(0, 1)
            j = random.randint(10, 99)
            k = random.randint(100, 199)
        else:
            i, j, k = None, None, None
        i = cOmm.bcast(i, root=mAster_rank)
        j = cOmm.bcast(j, root=mAster_rank)
        k = cOmm.bcast(k, root=mAster_rank)
        l = random.randint(2, 4)
        empty = sorted(random.sample(range(0, j), int(j / l)))

        sa = random.randint(1, 100)/300
        A = spspa.random(j, j, sa, format='lil')
        if i == 0:
            A[empty, :] = 0
        elif i == 1:
            A[:, empty] = 0
        else:
            raise Exception()
        if i == 0:
            A = A.tocsr()
        elif i == 1:
            A= A.tocsc()
        else:
            raise Exception()
        A0 = cOmm.gather(A, root=mAster_rank)
        if rAnk == mAster_rank:
            A0 = np.sum(A0)

        A = GlobalMatrix(A)

        A = A.___PRIVATE_gather_M_to_core___(clean_local=True, splitting_factor=k)
        if rAnk == mAster_rank:
            np.testing.assert_array_almost_equal(A.toarray(), A0.toarray())

    return 1

def test_TOOLS_NO7_linear_algebra_EWC_test():
    if rAnk == mAster_rank:
        print("::: [test_TOOLS_NO7_linear_algebra_EWC_test] ....", flush=True)

    for _ in range(2):
        if rAnk == mAster_rank:
            c = random.uniform(0, 0.15)
            i = random.randint(2, 3)
            j = random.randint(1, 3)
            k = random.randint(2, 4)
            l = random.randint(2, 3)
            m = random.randint(1, 3)
            n = random.randint(2, 4)
        else:
            c, i, j, k, l, m, n = None, None, None, None, None, None, None
        c, i, j, k, l, m, n = cOmm.bcast([c, i, j, k, l, m, n], root=mAster_rank)
        if c < 0.1: c = 0
        mesh = MeshGenerator('crazy', c=c)([f'Lobatto:{i}', f'Lobatto:{j}', f'Lobatto:{k}'], EDM='debug')
        space = SpaceInvoker('polynomials')([('Lobatto', l), ('Lobatto', m), ('Lobatto', n)])
        FC = FormCaller(mesh, space)

        f0 = FC('0-f', is_hybrid=False)
        f1 = FC('1-f', is_hybrid=False)
        f2 = FC('2-f', is_hybrid=False)

        E10 = f0.matrices.incidence
        E21 = f1.matrices.incidence
        E32 = f2.matrices.incidence


        E210 = E21 @ E10
        E321 = E32 @ E21

        for i in mesh.elements:
            e210 = E210[i].toarray()
            e321 = E321[i].toarray()
            np.testing.assert_almost_equal(np.sum(np.abs(e210)), 0)
            np.testing.assert_almost_equal(np.sum(np.abs(e321)), 0)

        M1 = f1.matrices.mass
        M2 = f2.matrices.mass

        blocks = ([M2  , E21],
                  [None, M1 ])

        BMAT = mif.bmat(blocks)

        for i in BMAT:
            bi = BMAT[i]
            m2 = M2[i]
            e21 = E21[i]
            m1 = M1[i]
            Bi = spspa.bmat(([m2, e21],[None, m1]), format='csc')
            np.testing.assert_almost_equal(np.sum(np.abs(bi - Bi)), 0)

    return 1

def test_TOOLS_NO8_GlobalMatrix_dot_product_test():
    """"""
    if rAnk == mAster_rank:
        print("... [test_TOOLS_NO8_GlobalMatrix_dot_product_test] ....", flush=True)

    for _ in range(20): # do multiple random tests.
        # A @ B @ C, A: (i,j), B: (j, k), C:(k, l)
        if rAnk == mAster_rank:
            i = random.randint(1, sIze*2)
            j = random.randint(1, sIze*2)
            k = random.randint(1, sIze*2)
            l = random.randint(1, sIze*2)

            type_A = random.randint(0, 1)
            type_B = random.randint(0, 1)
            type_C = random.randint(0, 1)
        else:
            i, j, k, l = None, None, None, None
            type_A, type_B, type_C = None, None, None
        i, j, k, l = cOmm.bcast([i, j, k, l], root=mAster_rank)
        type_A, type_B, type_C = cOmm.bcast([type_A, type_B, type_C], root=mAster_rank)
        type_A = ['row', 'col'][type_A]
        type_B = ['row', 'col'][type_B]
        type_C = ['row', 'col'][type_C]

        if type_A == 'row':
            A = ___generate_random_row_major_GM___(i, j)
        else:
            A = ___generate_random_col_major_GM___(i, j)

        if type_B == 'row':
            B = ___generate_random_row_major_GM___(j, k)
        else:
            B = ___generate_random_col_major_GM___(j, k)

        if type_C == 'row':
            C = ___generate_random_row_major_GM___(k, l)
        else:
            C = ___generate_random_col_major_GM___(k, l)

        ABC = A @ B @ C
        a = A.___PRIVATE_gather_M_to_core___()
        b = B.___PRIVATE_gather_M_to_core___()
        c = C.___PRIVATE_gather_M_to_core___()
        abc = ABC.___PRIVATE_gather_M_to_core___()

        if rAnk == mAster_rank:
            ___ = a @ b @ c
            np.testing.assert_almost_equal(np.sum(np.abs(abc - ___)), 0)

    return 1
def ___generate_random_row_major_GM___(i, j, s=None):
    """Make a random row major sparse matrix of shape (i,j) at sparsity=s.

    :param i:
    :param j:
    :param s:
    :return:
    """
    if s is None:
        s = random.uniform(0, 0.1)
        if s < 0.02: s = 0
    if rAnk == mAster_rank:

        random_list = random.sample(range(0, i), i)
        distribution = [i // sIze + (1 if x < i % sIze else 0) for x in range(sIze)]
        no_empty_rows = list()
        _ = 0
        for r in range(sIze):
            no_empty_rows.append(random_list[_:_+distribution[r]])
            _ += distribution[r]

    else:
        no_empty_rows = None

    no_empty_rows = cOmm.scatter(no_empty_rows, root=mAster_rank)

    _ = spspa.random(i, j, s, format='csr')

    A = spspa.lil_matrix((i,j))
    A[no_empty_rows,:] = _[no_empty_rows,:]
    A = A.tocsr()
    A = GlobalMatrix(A)
    A.IS.regularly_distributed = 'row'
    A.___PRIVATE_self_regularity_checker___()

    return A
def ___generate_random_col_major_GM___(i, j, s=None):
    """Make a random col major sparse matrix of shape (i,j) at sparsity=s.

    :param i:
    :param j:
    :param s:
    :return:
    """
    if s is None:
        s = random.uniform(0, 0.1)
        if s < 0.02: s = 0
    if rAnk == mAster_rank:

        random_list = random.sample(range(0, j), j)
        distribution = [j // sIze + (1 if x < j % sIze else 0) for x in range(sIze)]
        no_empty_cols = list()
        _ = 0
        for r in range(sIze):
            no_empty_cols.append(random_list[_:_+distribution[r]])
            _ += distribution[r]

    else:
        no_empty_cols = None

    no_empty_cols = cOmm.scatter(no_empty_cols, root=mAster_rank)

    _ = spspa.random(i, j, s, format='csc')

    A = spspa.lil_matrix((i,j))
    A[:,no_empty_cols] = _[:,no_empty_cols]
    A = A.tocsc()
    A = GlobalMatrix(A)
    A.IS.regularly_distributed = 'column'
    A.___PRIVATE_self_regularity_checker___()

    return A

def test_TOOLS_NO9_test_Chained_Gathering_Matrix():
    """"""
    if rAnk == mAster_rank:
        print("+++ [test_TOOLS_NO9_test_Chained_Gathering_Matrix] ....", flush=True)

    for _ in range(3):
        if rAnk == mAster_rank:
            i = random.randint(2, 3)
            j = random.randint(1, 3)
            k = random.randint(1, 2)
            l = random.randint(2, 3)
            m = random.randint(1, 2)
            n = random.randint(2, 3)

            mid = random.randint(0, 1)
        else:
            i, j, k, l, m, n = None, None, None, None, None, None
            mid = None
        i, j, k, l, m, n, mid = cOmm.bcast([i, j, k, l, m, n, mid], root=mAster_rank)

        MID = ('crazy', 'bridge_arch_cracked')[mid]

        mesh = MeshGenerator(MID)([f'Lobatto:{i}', f'Lobatto:{j}', f'Lobatto:{k}'], EDM='debug')
        space = SpaceInvoker('polynomials')([('Lobatto', l), ('Lobatto', m), ('Lobatto', n)])
        FC = FormCaller(mesh, space)
        f0 = FC('0-f', is_hybrid=False)
        f1 = FC('1-f', is_hybrid=False)
        f2 = FC('2-f', is_hybrid=False)
        f3 = FC('3-f', is_hybrid=False)

        GM0 = f0.numbering.gathering
        GM1 = f1.numbering.gathering
        GM2 = f2.numbering.gathering
        GM3 = f3.numbering.gathering

        CGM = Chain_Gathering_Matrix(GM0)
        for i in CGM:
            np.testing.assert_array_equal(CGM[i], GM0[i].full_vector)

        CGM0 = Chain_Gathering_Matrix([GM1, GM1])
        CGM1 = Chain_Gathering_Matrix([GM0, GM1, GM2])
        CGM2 = Chain_Gathering_Matrix([GM0, GM1, GM2])
        assert CGM1 == CGM2 and CGM1 is not CGM2, "Check equal!"
        CGM2 = Chain_Gathering_Matrix([GM1, GM2, GM3])
        assert CGM1 != CGM2
        assert CGM2 != CGM

        NUM0 = GM0.GLOBAL_num_dofs
        NUM1 = GM1.GLOBAL_num_dofs
        for i in CGM1:
            cgm = CGM1[i]
            gm0 = GM0[i].full_vector
            gm1 = GM1[i].full_vector
            gm2 = GM2[i].full_vector
            c_g_m = np.concatenate([gm0, gm1+NUM0, gm2+NUM0+NUM1])
            np.testing.assert_array_equal(cgm, c_g_m)

        assert CGM2.GLOBAL_num_dofs == GM1.GLOBAL_num_dofs + GM2.GLOBAL_num_dofs + GM3.GLOBAL_num_dofs

        CGM3 = Chain_Gathering_Matrix([GM0, GM1, GM2, GM3])
        CGM4 = Chain_Gathering_Matrix([GM0, GM0, GM2, GM3])
        CGM5 = Chain_Gathering_Matrix([GM0, GM3])
        CGM6 = Chain_Gathering_Matrix([GM2,])

        assert CGM0.mesh_type == '_3dCSCG'
        assert CGM1.mesh_type == '_3dCSCG'
        assert CGM2.mesh_type == '_3dCSCG'
        assert CGM3.mesh_type == '_3dCSCG'
        assert CGM4.mesh_type == '_3dCSCG'
        assert CGM5.mesh_type == '_3dCSCG'
        assert CGM6.mesh_type == '_3dCSCG'

        assert CGM6.GLOBAL_num_dofs == GM2.GLOBAL_num_dofs
        assert CGM0.GLOBAL_num_dofs == 2 * GM1.GLOBAL_num_dofs
        assert CGM4.GLOBAL_num_dofs == 2 * GM0.GLOBAL_num_dofs + GM2.GLOBAL_num_dofs + GM3.GLOBAL_num_dofs

        if rAnk == mAster_rank:
            COMM_CHECK = random.randint(3, 5)
            i = random.randint(0, 6)
        else:
            COMM_CHECK = None
            i = None
        i, COMM_CHECK = cOmm.bcast([i, COMM_CHECK], root=mAster_rank)
        CGM = (CGM0, CGM1, CGM2, CGM3, CGM4, CGM5, CGM6)[i]
        assert CGM.mesh_type == '_3dCSCG'

        gnd_S = list()
        for gm in CGM.GMs:
            gnd_S.append(gm.GLOBAL_num_dofs)

        for i in CGM:
            cgm = CGM[i]
            assert len(set(cgm)) == len(cgm), "Having repeated numbering? BAD BAD BAD!"
            gm_list = list()
            for j, gm in enumerate(CGM.GMs):
                gm_list.append(gm[i].full_vector+sum(gnd_S[0:j]))
            c_g_m = np.concatenate(gm_list)
            np.testing.assert_array_equal(cgm, c_g_m)


        GND = CGM.GLOBAL_num_dofs

        n = 0

        for m in range(GND):
            I = CGM.do.find.elements_contain_dof_numbered(m)
            if I is None:
                exclude_list = list()
            else:
                I = I[0]
                assert isinstance(I, list)
                exclude_list = I
                for i in I:
                    assert m in CGM[i], f"m not in element {i}!"

            for j in CGM:
                if j not in exclude_list:
                    assert m not in CGM[j], "m also in some elements else!"

            OUTPUT = CGM.do.find.elements_and_local_indices_of_dof(m)
            if OUTPUT is None:
                GET_it = 0
            else:
                GET_it = 1
                elements, indices = OUTPUT
                for e, ind in zip(elements, indices):
                    assert CGM[e][ind] == m, f"CGM[{e}][{ind}] = {CGM[e][ind]} which != {m}!"

            if n < COMM_CHECK and m % COMM_CHECK == 0:
                n += 1
                GET = cOmm.gather(GET_it, root=mAster_rank)
                if rAnk == mAster_rank:
                    assert sum(GET) >= 1, f"At least in one core we will find element(s) contain dof {m}!"

    return 1

def test_TOOLS_NO10_test_EWC_SparseMatrix_Customize():
    """"""
    if rAnk == mAster_rank:
        print("||| [test_TOOLS_NO10_test_EWC_SparseMatrix_Customize] ....", flush=True)

    for _ in range(1):
        if rAnk == mAster_rank:
            i = random.randint(2, 3)
            j = random.randint(1, 3)
            k = random.randint(1, 2)
            l = random.randint(2, 3)
            m = random.randint(1, 2)
            n = random.randint(2, 3)

            mid = random.randint(0, 1)
            ish = random.randint(0, 1)
        else:
            i, j, k, l, m, n = None, None, None, None, None, None
            mid = None
            ish = None
        i, j, k, l, m, n, mid, ish = cOmm.bcast([i, j, k, l, m, n, mid, ish], root=mAster_rank)

        MID = ('crazy', 'bridge_arch_cracked')[mid]
        ISH = (True, False)[ish]

        mesh = MeshGenerator(MID)([f'Lobatto:{i}', f'Lobatto:{j}', f'Lobatto:{k}'], EDM='debug')
        space = SpaceInvoker('polynomials')([('Lobatto', l), ('Lobatto', m), ('Lobatto', n)])
        FC = FormCaller(mesh, space)

        f2 = FC('2-f', is_hybrid=ISH)
        f3 = FC('3-f', is_hybrid=ISH)

        M3 = f3.matrices.mass
        E32 = f2.matrices.incidence

        f2_GND = f2.num.GLOBAL_dofs
        f3_GND = f3.num.GLOBAL_dofs
        GND = f2_GND + f3_GND
        BMAT = [[M3, E32],[E32.T, None]]
        SYSTEM = bmat(BMAT)
        SYSTEM.gathering_matrices=([f3,f2], [f3,f2])



        if rAnk == mAster_rank:
            AAA =  int(GND / 30)
            if AAA > 10:
                HOW_MANY_LINES = 10
            else:
                HOW_MANY_LINES = AAA
            III = random.sample(range(0, GND), HOW_MANY_LINES)

        else:
            III = None
        III = cOmm.bcast(III, root=mAster_rank)
        for iii in III:
            SYSTEM.customize.clear_global_row(iii)
        SYSTEM_ASSEMBLED = SYSTEM.assembled

        for jjj in III:
            assert SYSTEM_ASSEMBLED.M[jjj].nnz == 0
        M = SYSTEM_ASSEMBLED.___PRIVATE_gather_M_to_core___()
        for jjj in range(GND):
            if jjj not in III:
                if rAnk == mAster_rank:
                    assert M[jjj].nnz != 0

        SYSTEM = M3
        GND = f3_GND
        if rAnk == mAster_rank:
            AAA =  int(GND / 20)
            if AAA > 10:
                HOW_MANY_LINES = 10
            else:
                HOW_MANY_LINES = AAA
            III = random.sample(range(0, GND), HOW_MANY_LINES)
        else:
            III = None

        III = cOmm.bcast(III, root=mAster_rank)
        for iii in III:
            SYSTEM.customize.clear_global_row(iii)
        SYSTEM_ASSEMBLED = SYSTEM.assembled

        for jjj in III:
            assert SYSTEM_ASSEMBLED.M[jjj].nnz == 0
        M = SYSTEM_ASSEMBLED.___PRIVATE_gather_M_to_core___()
        for jjj in range(GND):
            if jjj not in III:
                if rAnk == mAster_rank:
                    assert M[jjj].nnz != 0

        f0 = FC('0-f', is_hybrid=False)
        f1 = FC('1-f', is_hybrid=False)

        M1 = f1.matrices.mass
        E10 = f0.matrices.incidence

        f0_GND = f0.num.GLOBAL_dofs
        f1_GND = f1.num.GLOBAL_dofs
        GND = f0_GND + f1_GND
        BMAT = [[M1, E10],[E10.T, None]]
        SYSTEM = bmat(BMAT)
        SYSTEM.gathering_matrices=([f1,f0], [f1,f0])

        if rAnk == mAster_rank:
            AAA =  int(GND / 30)
            if AAA > 10:
                HOW_MANY_LINES = 10
            else:
                HOW_MANY_LINES = AAA
            III = random.sample(range(0, GND), HOW_MANY_LINES)
        else:
            III = None
        III = cOmm.bcast(III, root=mAster_rank)
        for iii in III:
            SYSTEM.customize.clear_global_row(iii)
        SYSTEM_ASSEMBLED = SYSTEM.assembled

        for jjj in III:
            assert SYSTEM_ASSEMBLED.M[jjj].nnz == 0
        M = SYSTEM_ASSEMBLED.___PRIVATE_gather_M_to_core___()
        for jjj in range(GND):
            if jjj not in III:
                if rAnk == mAster_rank:
                    assert M[jjj].nnz != 0

        SYSTEM = bmat(BMAT)
        SYSTEM.gathering_matrices=([f1,f0], [f1,f0])
        if rAnk == mAster_rank:
            AAA =  int(GND / 30)
            if AAA > 10:
                HOW_MANY_LINES = 10
            else:
                HOW_MANY_LINES = AAA
            III = random.sample(range(0, GND), HOW_MANY_LINES)
            VVV = np.random.rand(HOW_MANY_LINES)
        else:
            III, VVV = None, None
        III, VVV = cOmm.bcast([III, VVV], root=mAster_rank) # we will set the value at M[i,j], i in III, j in JJJ.

        for i, v in zip(III, VVV):
            SYSTEM.customize.set_assembled_M_ij_to(i, i, v)
        SYSTEM_ASSEMBLED = SYSTEM.assembled
        M = SYSTEM_ASSEMBLED.___PRIVATE_gather_M_to_core___()

        if rAnk == mAster_rank:
            for i, v in zip(III, VVV):
                assert M[i,i] == v

        SYSTEM = bmat(BMAT)
        SYSTEM.gathering_matrices=([f1,f0], [f1,f0])
        if rAnk == mAster_rank:
            AAA =  int(GND / 30)
            if AAA > 10:
                HOW_MANY_LINES = 10
            else:
                HOW_MANY_LINES = AAA
            III = random.sample(range(0, GND), HOW_MANY_LINES)
        else:
            III, VVV = None, None
        III = cOmm.bcast(III, root=mAster_rank) # we will set the value at M[i,j], i in III, j in JJJ.

        for i in III:
            SYSTEM.customize.identify_global_row(i)
        SYSTEM_ASSEMBLED = SYSTEM.assembled
        M = SYSTEM_ASSEMBLED.___PRIVATE_gather_M_to_core___()

        if rAnk == mAster_rank:
            for i in III:
                assert M[i].nnz == 1
                assert M[i,i] == 1


    return 1






def ___runner_test_function_11___(t0, dt, steps, A=10):
    """
    Parameters
    ----------
    t0 :
    dt :
    steps :
    A : A key ward input.

    Returns
    -------
    o1 :
    o2 :
    """
    SI = SimpleIterator(t0=t0, dt=dt, max_steps=steps)
    SI(___TEST_SOLVER___, [0, 0])
    SI.run()
    if rAnk == mAster_rank:
        result = np.array(SI.RDF.to_numpy()[-1,:][-2:], dtype=float)
        return result[0]+A, result[1]-A

def test_TOOLS_NO11_test_ParallelMatrix3dInputRunner():
    """"""
    if rAnk == mAster_rank:
        print("->- [test_TOOLS_NO11_test_ParallelMatrix3dInputRunner] ....", flush=True)

        i1 = random.randint(2,3)
        i2 = random.randint(1,4)
        J = random.randint(1,3)

        T0 = (random.sample(range(0, 2*i1), i1),
              random.sample(range(2*i1+1, 2*i1+1+2*i2), i2))
        DT = ([(random.random()+0.1) * np.pi /10 for _ in range(i1)],
              [(random.random()+0.1) * np.sqrt(2) /10 for _ in range(i2)])
        STEPS = random.sample(range(1, 2*J), J)

        A = random.uniform(0,10)
    else:
        T0, DT, STEPS, A = None, None, None, None

    T0, DT, STEPS, A = cOmm.bcast([T0, DT, STEPS, A], root=mAster_rank)


    PR1 = ParallelMatrix3dInputRunner(___runner_test_function_11___)
    PR1.iterate(T0, DT, STEPS, writeto='pmr_test.txt', A=A)

    if rAnk == mAster_rank:

        DR1 = len(PR1._SR_.rdf)
        if DR1 < 5: # if we have less than 5 runs, we re-run all.
            DEL_ROWS = DR1
        elif DR1 > 15:
            DEL_ROWS = 0 # re-run none.
        elif DR1 > 11:
            DEL_ROWS = random.randint(0, int(DR1/3))
        else:
            DEL_ROWS = random.randint(0, DR1)

        with open('pmr_test.txt', 'r') as f:
            contents = f.readlines()
            LEN1 = len(contents)
            ROWS = random.sample(range(0, LEN1-1), random.randint(1,10))
            ROWS.extend([-1,-2,-3])
            TO_BE_CHECK = list()
            for r in ROWS:
                TO_BE_CHECK.append(contents[r])

        if DEL_ROWS != 0:
            with open('pmr_test.txt', 'w') as f:
                    for con in contents[:-DEL_ROWS]:
                        f.write(con)

    PR2 = ParallelMatrix3dInputRunner(___runner_test_function_11___)
    PR2.iterate(T0, DT, STEPS, writeto='pmr_test.txt', A=A)

    if rAnk == mAster_rank:
        with open('pmr_test.txt', 'r') as f:
            contents = f.readlines()
            for i, r in enumerate(ROWS):

                if len(TO_BE_CHECK[i]) >= 138:
                    assert TO_BE_CHECK[i][:138] == contents[r][:138]
                else:
                    assert TO_BE_CHECK[i] == contents[r]

    PR3 = RunnerDataReader('pmr_test.txt')


    if rAnk == mAster_rank:
        D = PR2._SR_.rdf.to_numpy()
    else:
        D = None

    D = cOmm.bcast(D, root=mAster_rank)
    assert np.all(D == PR3.results.to_numpy())
    assert PR3.___lock_iterate___

    if rAnk == mAster_rank:
        os.remove('pmr_test.txt')

    return 1



def test_TOOLS_NO12_EWC_assembling_test():
    """"""
    if rAnk == mAster_rank:
        print("AAA [test_TOOLS_NO12_EWC_assembling_test] ....", flush=True)

    if rAnk == mAster_rank:
        el1 = random.randint(1,3)
        el2 = random.randint(1,3)
        el3 = random.randint(1,3)
        c = random.uniform(0.0, 0.3)
        if c < 0.1:c = 0
    else:
        el1, el2, el3, c = [None for _ in range(4)]
    el1, el2, el3, c = cOmm.bcast([el1, el2, el3, c], root=mAster_rank)

    mesh = MeshGenerator('crazy', c=c)([el1, el2, el3])
    space = SpaceInvoker('polynomials')([('Lobatto', el3), ('Lobatto', el2), ('Lobatto', el1)])
    FC = FormCaller(mesh, space)
    f0 = FC('0-f', is_hybrid=False)
    f1 = FC('1-f', is_hybrid=False)
    f2 = FC('2-f', is_hybrid=False)
    f3 = FC('2-f', is_hybrid=False)

    M0 = f0.matrices.mass
    M1 = f1.matrices.mass
    M2 = f2.matrices.mass
    M3 = f3.matrices.mass
    GM0 = f0.numbering.gathering
    GM1 = f1.numbering.gathering
    GM2 = f2.numbering.gathering
    GM3 = f3.numbering.gathering

    aM0 = ___brutal_force_EWC_matrix_assembling___(M0, GM0, GM0)
    aM1 = ___brutal_force_EWC_matrix_assembling___(M1, GM1, GM1)
    aM2 = ___brutal_force_EWC_matrix_assembling___(M2, GM2, GM2)
    aM3 = ___brutal_force_EWC_matrix_assembling___(M3, GM3, GM3)
    M0.gathering_matrices = (GM0, GM0)
    M1.gathering_matrices = (GM1, GM1)
    M2.gathering_matrices = (GM2, GM2)
    M3.gathering_matrices = (GM3, GM3)
    AM0 = M0.assembled
    AM1 = M1.assembled
    AM2 = M2.assembled
    AM3 = M3.assembled
    np.testing.assert_array_almost_equal(aM0.toarray(), AM0.M.toarray())
    np.testing.assert_array_almost_equal(aM1.toarray(), AM1.M.toarray())
    np.testing.assert_array_almost_equal(aM2.toarray(), AM2.M.toarray())
    np.testing.assert_array_almost_equal(aM3.toarray(), AM3.M.toarray())

    def p(t, x, y, z): return -0.001 + x + 2.1254*y + 0.1234*z + t
    scalar = FC('scalar', p)

    f0.TW.func.do.set_func_body_as(scalar)
    f0.TW.do.push_all_to_instant(0)
    f0.discretize()
    c0 = f0.cochain.local
    aC0 = ___brutal_force_EWC_vector_assembling___(c0, GM0)
    C0 = f0.cochain.EWC
    AC0 = C0.assembled
    np.testing.assert_array_almost_equal(aC0.toarray(), AC0.V.toarray())

    mesh = MeshGenerator2D('crazy', c=c)([el1, el2])
    space = SpaceInvoker2D('polynomials')([('Lobatto', el3), ('Lobatto', el2)])
    FC = FormCaller2D(mesh, space)
    f0 = FC('0-f-i', is_hybrid=False)
    f1 = FC('1-f-i', is_hybrid=False)
    f2 = FC('2-f-i', is_hybrid=False)
    M0 = f0.matrices.mass
    M1 = f1.matrices.mass
    M2 = f2.matrices.mass
    GM0 = f0.numbering.gathering
    GM1 = f1.numbering.gathering
    GM2 = f2.numbering.gathering
    aM0 = ___brutal_force_EWC_matrix_assembling___(M0, GM0, GM0)
    aM1 = ___brutal_force_EWC_matrix_assembling___(M1, GM1, GM1)
    aM2 = ___brutal_force_EWC_matrix_assembling___(M2, GM2, GM2)
    M0.gathering_matrices = (GM0, GM0)
    M1.gathering_matrices = (GM1, GM1)
    M2.gathering_matrices = (GM2, GM2)
    AM0 = M0.assembled
    AM1 = M1.assembled
    AM2 = M2.assembled
    np.testing.assert_array_almost_equal(aM0.toarray(), AM0.M.toarray())
    np.testing.assert_array_almost_equal(aM1.toarray(), AM1.M.toarray())
    np.testing.assert_array_almost_equal(aM2.toarray(), AM2.M.toarray())

    E10 = f0.matrices.incidence
    E10.gathering_matrices = (GM1, GM0)
    M = bmat(([M0],[E10]))

    aM = ___brutal_force_EWC_matrix_assembling___(M, *M.gathering_matrices)
    AM = M.assembled
    np.testing.assert_array_almost_equal(aM.toarray(), AM.M.toarray())

    def p(t, x, y): return -0.001 + x + 2.1254*y + t
    scalar = FC('scalar', p)

    f0.TW.func.do.set_func_body_as(scalar)
    f0.TW.do.push_all_to_instant(0)
    f0.discretize()
    c0 = f0.cochain.local
    aC0 = ___brutal_force_EWC_vector_assembling___(c0, GM0)
    C0 = f0.cochain.EWC
    AC0 = C0.assembled
    np.testing.assert_array_almost_equal(aC0.toarray(), AC0.V.toarray())

    return 1

def ___brutal_force_EWC_matrix_assembling___(EWC, GM0, GM1):
    """Do the brutal force assembling."""
    N0 = GM0.GLOBAL_num_dofs
    N1 = GM1.GLOBAL_num_dofs

    AM = spspa.lil_matrix((N0, N1))
    for i in EWC:
        Mi = EWC[i]
        gm0 = GM0[i]
        gm1 = GM1[i]

        for j,m in enumerate(gm0):
            for k,n in enumerate(gm1):
                AM[m,n] += Mi[j,k]

    return AM
def ___brutal_force_EWC_vector_assembling___(EWC, GM):
    """Do the brutal force assembling."""
    N = GM.GLOBAL_num_dofs

    A = spspa.lil_matrix((N, 1))
    for i in EWC:
        Vi = EWC[i]
        gm = GM[i]

        for j,m in enumerate(gm):
            A[m,0] += Vi[j]

    return A





def test_TOOLS_NO13_EWC_Customize_CSCG_partial_dofs():
    """"""
    if rAnk == mAster_rank:
        print("-p- [test_TOOLS_NO13_EWC_Customize_CSCG_partial_dofs] ....", flush=True)

    if rAnk == mAster_rank:
        ISH = (True, False, False)[random.randint(0, 2)]
        LOAD = random.randint(25, 75)
    else:
        ISH = None
        LOAD = None
    ISH, LOAD = cOmm.bcast([ISH, LOAD], root=mAster_rank)

    mesh, space = random_mesh_and_space_of_total_load_around(LOAD,
             exclude_periodic=True,
             domain_boundary_distribution_regularities='Regular:interfaces-not-shared-by-regions')

    #-------- test partial_dofs only having boundary dofs involved====
    bns = mesh.boundaries.names
    if rAnk == mAster_rank:
        BNS = list()
        for _ in range(4):
            i = random.randint(1, len(bns)) # how many boundaries to be included.
            BNS.append(random.sample(bns, i))
    else:
        BNS = None
    BNS = cOmm.bcast(BNS, root=mAster_rank)

    FC = FormCaller(mesh, space)
    f0 = FC('0-f', is_hybrid=ISH)
    f1 = FC('1-f', is_hybrid=ISH)

    pd0 = PartialDofs(f0)
    pd0.include.boundaries(BNS[0])

    M1 = f1.matrices.mass
    E10 = f0.matrices.incidence
    BMAT = [[M1, E10], [E10.T, None]]
    SYSTEM = bmat(BMAT)
    SYSTEM.gathering_matrices=([f1,f0], [f1,f0])

    SYSTEM.customize.identify_global_rows_according_to_CSCG_partial_dofs(1, pd0)

    SYSTEM = SYSTEM.assembled
    M = SYSTEM.___PRIVATE_gather_M_to_core___()

    Gbd = f0.numbering.GLOBAL_boundary_dofs

    if rAnk == mAster_rank:
        bd = list()
        for bn in BNS[0]:
            bd.extend(Gbd[bn])

        bd = np.array(bd) + f1.num.GLOBAL_dofs

        for i in range(f1.num.GLOBAL_dofs + f0.num.GLOBAL_dofs):
            if i in bd:
                assert M[i].nnz == 1, f"M[i].nnz={M[i].nnz}"
            else:
                assert M[i].nnz > 1

    # ------------ test for f2 -------------------------------------------------
    mesh, space = random_mesh_and_space_of_total_load_around(int(LOAD * 2 / 3),
                                                             exclude_periodic=True)

    #-------- test partial_dofs only having boundary dofs involved===
    bns = mesh.boundaries.names
    if rAnk == mAster_rank:
        BNS = list()
        for _ in range(4):
            i = random.randint(1, len(bns)) # how many boundaries to be included.
            BNS.append(random.sample(bns, i))
    else:
        BNS = None
    BNS = cOmm.bcast(BNS, root=mAster_rank)

    FC = FormCaller(mesh, space)
    f2 = FC('2-f', is_hybrid=ISH)
    f3 = FC('3-f', is_hybrid=ISH)
    pd2 = PartialDofs(f2)
    pd2.include.boundaries(BNS[2])

    M3 = f3.matrices.mass
    E32 = f2.matrices.incidence
    BMAT = [[M3, E32], [E32.T, None]]
    SYSTEM = bmat(BMAT)
    SYSTEM.gathering_matrices=([f3,f2], [f3,f2])
    SYSTEM.customize.identify_global_rows_according_to_CSCG_partial_dofs(1, pd2)
    SYSTEM = SYSTEM.assembled
    M = SYSTEM.___PRIVATE_gather_M_to_core___()

    Gbd_f2 = f2.numbering.GLOBAL_boundary_dofs

    if rAnk == mAster_rank:
        bdf2 = list()
        for bn in BNS[2]:
            bdf2.extend(Gbd_f2[bn])

        bdf2 = np.array(bdf2) + f3.num.GLOBAL_dofs

        for i in bdf2:
            assert M[i].nnz == 1, f"M[i].nnz={M[i].nnz}"
            assert M[i,i] == 1

    t2 = FC('2-t')
    pdt = PartialDofs(t2)
    pdt.include.boundaries(BNS[2])

    M2 = f2.matrices.mass
    E32 = f2.matrices.incidence
    T2 = t2.matrices.trace
    BMAT = [[M2, E32.T, T2.T],
            [E32, None, None],
            [T2, None, None]]
    SYSTEM = bmat(BMAT)
    SYSTEM.gathering_matrices = ([f2, f3, t2], [f2, f3, t2])

    SYSTEM.customize.off_diagonally_identify_rows_according_to_two_CSCG_partial_dofs(
        2, 0, pdt, pd2)
    SYSTEM = SYSTEM.assembled
    M = SYSTEM.___PRIVATE_gather_M_to_core___()


    GM_t2 = t2.numbering.gathering
    GM_F2 = f2.numbering.gathering

    gn_T2 = dict()
    gn_F2 = dict()
    for e in pdt:
        gn_T2[e] = GM_t2[e][pdt.interpreted_as.local_dofs[e]]
        gn_F2[e] = GM_F2[e][pd2.interpreted_as.local_dofs[e]]


    gn_T2 = cOmm.gather(gn_T2, root=mAster_rank)
    gn_F2 = cOmm.gather(gn_F2, root=mAster_rank)

    if rAnk == mAster_rank:
        AT2 = dict()
        for _ in gn_T2: AT2.update(_)
        AF2 = dict()
        for _ in gn_F2: AF2.update(_)

        for e in AT2:
            gt = AT2[e]
            gf = AF2[e]

            gt = gt + f2.num.GLOBAL_dofs + f3.num.GLOBAL_dofs

            for i, j in zip(gt, gf):
                assert M[i].nnz == 1 and M[i,j] == 1

    # --------- single EWC -----------------------------------------------------------
    BMAT = [[M2,],]
    SYSTEM = bmat(BMAT)
    SYSTEM.gathering_matrices = ([f2, ], [f2, ])
    SYSTEM.customize.off_diagonally_identify_rows_according_to_two_CSCG_partial_dofs(
        0, 0, pd2, pd2)
    SYSTEM = SYSTEM.assembled
    M = SYSTEM.___PRIVATE_gather_M_to_core___()

    if rAnk == mAster_rank:
        # noinspection PyUnboundLocalVariable
        for e in AF2:
            gf = AF2[e]
            for i in gf:
                assert M[i].nnz == 1 and M[i,i] == 1
    # below, we check M2 is not changed by above adjusting.
    pd2_new = PartialDofs(f2)
    pd2_new.include.boundaries(BNS[3])
    SYSTEM = bmat(BMAT)
    SYSTEM.gathering_matrices = ([f2, ], [f2, ])
    SYSTEM.customize.off_diagonally_identify_rows_according_to_two_CSCG_partial_dofs(
        0, 0, pd2_new, pd2_new)
    SYSTEM = SYSTEM.assembled
    M_new = SYSTEM.___PRIVATE_gather_M_to_core___()
    gn_F2 = dict()
    for e in pd2_new:
        gn_F2[e] = GM_F2[e][pd2_new.interpreted_as.local_dofs[e]]
    gn_F2 = cOmm.gather(gn_F2, root=mAster_rank)
    if rAnk == mAster_rank:
        AF2_new = dict()
        for _ in gn_F2: AF2_new.update(_)
        ALL = list()
        for e in AF2_new:
            ALL.extend(AF2_new[e])
        for i in range(np.shape(M_new)[0]):
            if i in ALL:
                assert M_new[i].nnz == 1 and M_new[i,i] == 1
            else:
                assert M_new[i].nnz > 1

    return 1


def test_TOOLS_NO14_partial_cochain_with_3dCSCG_form_BC():
    """"""
    if rAnk == mAster_rank:
        print("-C- [test_TOOLS_NO14_partial_cochain_with_3dCSCG_form_BC] ....", flush=True)

    def Pressure(t, x, y, z): return 2.1 + t + np.cos(np.pi * x) * np.cos(2 * np.pi * y) * np.cos(3 * np.pi * z)
    def velocity_x(t, x, y, z): return 2.1 + t + np.cos(1.5*np.pi * x) * np.cos(2.5 * np.pi * y) * np.cos(3.5 * np.pi * z)
    def velocity_y(t, x, y, z): return 2.1 + t + np.cos(0.5*np.pi * x) * np.cos(2.1 * np.pi * y) * np.cos(3 * np.pi * z)
    def velocity_z(t, x, y, z): return 2.1 + t + np.cos(0.8*np.pi * x) * np.cos(1.3 * np.pi * y) * np.cos(0.7 * np.pi * z)

    if rAnk == mAster_rank:
        ISH = (True, False, False, False)[random.randint(0, 3)]
        LOAD = random.randint(50, 85)
        time = random.random()
    else:
        ISH = None
        LOAD = None
        time = None
    ISH, LOAD, time = cOmm.bcast([ISH, LOAD, time], root=mAster_rank)

    mesh, space = random_mesh_and_space_of_total_load_around(LOAD, exclude_periodic=True)

    bns = mesh.boundaries.names
    if rAnk == mAster_rank:
        i = random.randint(1, len(bns)) # how many boundaries to be included.
        BNS = random.sample(bns, i)
    else:
        BNS = None
    BNS = cOmm.bcast(BNS, root=mAster_rank)
    bcDs = dict()
    for bn in BNS: bcDs[bn] = Pressure
    bcDv = dict()
    for bn in BNS: bcDv[bn] = [velocity_x, velocity_y, velocity_z]

    FC = FormCaller(mesh, space)

    BS = FC('scalar', bcDs)
    SS = FC('scalar', Pressure)

    BV = FC('vector', bcDv)
    SV = FC('vector', [velocity_x, velocity_y, velocity_z])

    #----  with 3d CSCG 0-form ---------------------------------------------------------------------
    f0 = FC('0-f', is_hybrid=ISH)
    f0.TW.BC.body = BS
    f0.TW.do.push_BC_to_instant(time)
    f0.BC.valid_boundaries = BNS
    f0pc = f0.BC.partial_cochain
    xi_et_sg = np.meshgrid(*space.nodes, indexing='ij')
    for i in f0pc:
        element = mesh.elements[i]
        local_dofs = f0pc.dofs.interpreted_as.local_dofs[i]
        local_cochain = f0pc.cochain[i]
        x, y, z = element.coordinate_transformation.mapping(*xi_et_sg)
        x = x.ravel('F')[local_dofs]
        y = y.ravel('F')[local_dofs]
        z = z.ravel('F')[local_dofs]
        v_exact = Pressure(time, x, y, z)
        np.testing.assert_array_almost_equal(local_cochain, v_exact)

    f0.TW.BC.body = SS
    f0.TW.do.push_BC_to_instant(time)
    f0.BC.valid_boundaries = BNS
    f0pc = f0.BC.partial_cochain
    xi_et_sg = np.meshgrid(*space.nodes, indexing='ij')
    for i in f0pc:
        element = mesh.elements[i]
        local_dofs = f0pc.dofs.interpreted_as.local_dofs[i]
        local_cochain = f0pc.cochain[i]
        x, y, z = element.coordinate_transformation.mapping(*xi_et_sg)
        x = x.ravel('F')[local_dofs]
        y = y.ravel('F')[local_dofs]
        z = z.ravel('F')[local_dofs]
        v_exact = Pressure(time, x, y, z)
        np.testing.assert_array_almost_equal(local_cochain, v_exact)

    #----  with 3d CSCG 2-form ---------------------------------------------------------------------
    f2 = FC('2-f', is_hybrid=ISH)
    f2.TW.func.body = SV
    f2.TW.do.push_func_to_instant(time)
    f2.discretize()
    f2_cochain = f2.cochain.local

    f2.TW.BC.body = BV
    f2.TW.do.push_BC_to_instant(time)
    f2.BC.valid_boundaries = BNS
    f2pc = f2.BC.partial_cochain
    for i in f2pc:
        local_dofs = f2pc.dofs.interpreted_as.local_dofs[i]
        local_cochain = f2pc.cochain[i]
        cochain_exact = f2_cochain[i][local_dofs]
        np.testing.assert_array_almost_equal(local_cochain, cochain_exact)

    f2.TW.BC.body = SV
    f2.TW.do.push_BC_to_instant(time)
    f2.BC.valid_boundaries = BNS
    f2pc = f2.BC.partial_cochain
    for i in f2pc:
        local_dofs = f2pc.dofs.interpreted_as.local_dofs[i]
        local_cochain = f2pc.cochain[i]
        cochain_exact = f2_cochain[i][local_dofs]
        np.testing.assert_array_almost_equal(local_cochain, cochain_exact)

    #------ with 3d CSCG 2-trace-form -----------------------------------------------
    t2 = FC('2-t')
    t2.TW.func.body = SV
    t2.TW.do.push_func_to_instant(time)
    t2.discretize()
    t2_cochain = t2.cochain.local
    # 1: standard vector
    t2.TW.BC.body = SV
    t2.TW.do.push_BC_to_instant(time)
    t2.BC.valid_boundaries = BNS
    t2pc = t2.BC.partial_cochain
    for i in t2pc:
        local_dofs = t2pc.dofs.interpreted_as.local_dofs[i]
        local_cochain = t2pc.cochain[i]
        cochain_exact = t2_cochain[i][local_dofs]
        np.testing.assert_array_almost_equal(local_cochain, cochain_exact)
    # 2: boundary-wise vector
    t2.TW.BC.body = BV
    t2.TW.do.push_BC_to_instant(time)
    t2.BC.valid_boundaries = BNS
    t2pc = t2.BC.partial_cochain
    for i in t2pc:
        local_dofs = t2pc.dofs.interpreted_as.local_dofs[i]
        local_cochain = t2pc.cochain[i]
        cochain_exact = t2_cochain[i][local_dofs]
        np.testing.assert_array_almost_equal(local_cochain, cochain_exact)

    t2.TW.func.body = SS
    t2.TW.do.push_func_to_instant(time)
    t2.discretize()
    t2_cochain = t2.cochain.local
    # 1: standard scalar
    t2.TW.BC.body = SS
    t2.TW.do.push_BC_to_instant(time)
    t2.BC.valid_boundaries = BNS
    t2pc = t2.BC.partial_cochain
    for i in t2pc:
        local_dofs = t2pc.dofs.interpreted_as.local_dofs[i]
        local_cochain = t2pc.cochain[i]
        cochain_exact = t2_cochain[i][local_dofs]
        np.testing.assert_array_almost_equal(local_cochain, cochain_exact)
    # 2: boundary-wise scalar
    t2.TW.BC.body = BS
    t2.TW.do.push_BC_to_instant(time)
    t2.BC.valid_boundaries = BNS
    t2pc = t2.BC.partial_cochain
    for i in t2pc:
        local_dofs = t2pc.dofs.interpreted_as.local_dofs[i]
        local_cochain = t2pc.cochain[i]
        cochain_exact = t2_cochain[i][local_dofs]
        np.testing.assert_array_almost_equal(local_cochain, cochain_exact)


    # test set_entries_according_to_CSCG_partial_cochains for EWC vectors.
    cf0 = EWC_ColumnVector(mesh, f0.num.basis)
    cf2 = EWC_ColumnVector(mesh, f2.num.basis)
    ct2 = EWC_ColumnVector(mesh, t2.num.basis)

    b = concatenate([cf0, cf2, ct2])
    b.gathering_matrix = [f0, f2, t2]
    b.customize.set_entries_according_to_CSCG_partial_cochains(2, t2pc)
    B = b.assembled
    B = B.___PRIVATE_gather_V_to_core___()

    GM = t2.numbering.gathering
    gn_t2 = dict()
    for e in t2pc:
        gn_t2[e] = GM[e][t2pc.dofs.interpreted_as.local_dofs[e]]
    gn_t2 = cOmm.gather(gn_t2, root=mAster_rank)
    if rAnk == mAster_rank:
        At2 = dict()
        for _ in gn_t2: At2.update(_)

        ALL = list()
        for e in At2:
            ALL.extend(At2[e])
        ALL = np.array(ALL) + f0.num.GLOBAL_dofs + f2.num.GLOBAL_dofs

        for i, bi in enumerate(B):
            if i in ALL:
                assert bi != 0
            else:
                assert bi == 0

    b.customize.set_constant_entries_according_to_CSCG_partial_dofs(2, t2pc, -1)
    B = b.assembled
    B = B.___PRIVATE_gather_V_to_core___()
    if rAnk == mAster_rank:
        for i, bi in enumerate(B):
            # noinspection PyUnboundLocalVariable
            if i in ALL:
                assert bi == -1
            else:
                assert bi == 0

    b = concatenate([ct2, cf2, cf0])
    b.gathering_matrix = [t2, f2, f0]
    b.customize.set_constant_entries_according_to_CSCG_partial_dofs(0, t2pc, -2.1415)
    B = b.assembled
    B = B.___PRIVATE_gather_V_to_core___()
    if rAnk == mAster_rank:
        ALL = ALL - f0.num.GLOBAL_dofs - f2.num.GLOBAL_dofs
        for i, bi in enumerate(B):
            if i in ALL:
                assert bi == -2.1415
            else:
                assert bi == 0



    return 1



def test_TOOLS_NO15_linear_system_apply_BC():
    """"""
    if rAnk == mAster_rank:
        LOAD = random.randint(50, 300)
        time = random.random()
        print(f"-S- [test_TOOLS_NO15_linear_system_apply_BC] @ load = {LOAD}, time=%.2f..."%time, flush=True)
    else:
        LOAD = None
        time = None
    LOAD, time = cOmm.bcast([LOAD, time], root=mAster_rank)

    mesh, space = random_mesh_and_space_of_total_load_around(LOAD, exclude_periodic=True, mesh_boundary_num='>=2')
    FC = FormCaller(mesh, space)

    def Pressure(t, x, y, z): return 2.5 + t + np.cos(np.pi * x) * np.cos(2 * np.pi * y) * np.cos(3 * np.pi * z)
    def velocity_x(t, x, y, z): return 2.5 + t + np.cos(1.5*np.pi * x) * np.cos(2.5 * np.pi * y) * np.cos(3.5 * np.pi * z)
    def velocity_y(t, x, y, z): return 2.5 + t + np.cos(0.5*np.pi * x) * np.cos(2.1 * np.pi * y) * np.cos(3.2 * np.pi * z)
    def velocity_z(t, x, y, z): return 2.5 + t + np.cos(0.8*np.pi * x) * np.cos(1.3 * np.pi * y) * np.cos(0.7 * np.pi * z)
    # SS = FC('scalar', Pressure)
    # SV = FC('vector', [velocity_x, velocity_y, velocity_z])

    bns = mesh.boundaries.names
    if rAnk == mAster_rank:
        i = random.randint(1, len(bns)-1) # how many boundaries to be included.
        BNS = random.sample(bns, i)
    else:
        BNS = None
    BNS = cOmm.bcast(BNS, root=mAster_rank)
    bcDs = dict()
    for bn in BNS: bcDs[bn] = Pressure
    # bcDv = dict()
    # for bn in BNS: bcDv[bn] = [velocity_x, velocity_y, velocity_z]
    BS = FC('scalar', bcDs)
    # BV = FC('vector', bcDv)


    BNS_com = list()
    for bn in bns:
        if bn not in BNS:
            BNS_com.append(bn)
    # bcDs_com = dict()
    # for bn in BNS_com: bcDs_com[bn] = Pressure
    bcDv_com = dict()
    for bn in BNS_com: bcDv_com[bn] = [velocity_x, velocity_y, velocity_z]
    # BS_com = FC('scalar', bcDs_com)
    BV_com = FC('vector', bcDv_com)


    # Poisson hybrid system ---------------------------------------------------------------------
    f2 = FC('2-f', is_hybrid=True)
    f3 = FC('3-f', is_hybrid=True)
    t2 = FC('2-t')

    M2 = f2.matrices.mass
    E32 = f2.matrices.incidence
    T2 = t2.matrices.trace
    BMAT = [[M2, E32.T, T2.T],
            [E32, None, None],
            [T2, None, None]]
    A = bmat(BMAT)
    A.gathering_matrices = ([f2, f3, t2], [f2, f3, t2])

    b0 = EWC_ColumnVector(mesh, f2.num.basis)
    b1 = EWC_ColumnVector(mesh, f3.num.basis)
    b2 = EWC_ColumnVector(mesh, t2.num.basis)
    b = concatenate([b0, b1, b2])
    b.gathering_matrix = [f2, f3, t2]

    Axb = LinearSystem(A, b)


    t2.TW.BC.body = BS
    t2.TW.do.push_BC_to_instant(time)
    t2.BC.valid_boundaries = BNS
    t2BC1 = t2.BC.partial_cochain

    Axb.customize.apply_strong_BC(2, 2, t2BC1)

    aA, ab = Axb.assembled
    aA = aA.___PRIVATE_gather_M_to_core___()
    ab = ab.___PRIVATE_gather_V_to_core___()


    GM_t2 = t2.numbering.gathering

    t2GBD = t2.numbering.GLOBAL_boundary_dofs
    if rAnk == mAster_rank:


        dofs_changed = list() # the rows that been changed in the global matrix.
        for bn in BNS:
            dofs_changed.extend(t2GBD[bn])
        dofs_changed = np.array(dofs_changed) + f2.num.GLOBAL_dofs \
                       + f3.num.GLOBAL_dofs
        for i in dofs_changed:
            assert aA[i].nnz == 1 and aA[i,i] == 1 and ab[i] != 0


        NOT_changed = list() # the rows that not been changed in the global matrix.
        for i in range(aA.shape[0]):
            if i in dofs_changed:
                pass
            else:
                NOT_changed.append(i)

        ab_SUB = ab[NOT_changed]
        assert np.all(ab_SUB==0), f"not change places must be all zero!"
        aA_SUB = aA[NOT_changed,:]
        aA_SUB = aA_SUB[f2.num.GLOBAL_dofs:, f2.num.GLOBAL_dofs:]
        assert aA_SUB.nnz == 0, f"not change places must be all zero!"


    gn_T2 = dict()
    for e in t2BC1:
        gn_T2[e] = GM_t2[e][t2BC1.dofs.interpreted_as.local_dofs[e]]

    gn_T2 = cOmm.gather(gn_T2, root=mAster_rank)

    if rAnk == mAster_rank:
        AT2 = dict()
        for _ in gn_T2: AT2.update(_)

        changed = list()

        for e in AT2:
            gt = AT2[e]
            gt = gt + f2.num.GLOBAL_dofs + f3.num.GLOBAL_dofs

            for i in gt:
                assert aA[i].nnz == 1 and aA[i,i] == 1 and ab[i] != 0

            changed.extend(gt)

        not_changed = list()
        for i in range(aA.shape[0]):
            if i in changed:
                pass
            else:
                not_changed.append(i)

        # noinspection PyUnboundLocalVariable
        assert not_changed == NOT_changed

        ab = ab[not_changed]
        assert np.all(ab==0), f"not change places must be all zero!"
        aA = aA[not_changed,:]
        aA = aA[f2.num.GLOBAL_dofs:, f2.num.GLOBAL_dofs:]
        assert aA.nnz == 0, f"not change places must be all zero!"

    t2BC2 = PartialDofs(t2)
    t2BC2.include.boundaries(BNS_com)
    f2.TW.BC.body = BV_com
    f2.TW.do.push_BC_to_instant(time)
    f2.BC.valid_boundaries = BNS_com
    f2BC2 = f2.BC.partial_cochain

    Axb.customize.apply_strong_BC(2, 0, t2BC2, f2BC2)
    aA, ab = Axb.assembled
    aA = aA.___PRIVATE_gather_M_to_core___()
    ab = ab.___PRIVATE_gather_V_to_core___()

    if rAnk == mAster_rank:


        dofs_changed = list() # the rows that been changed in the global matrix.
        for bn in t2GBD: # all the boundary dofs.
            dofs_changed.extend(t2GBD[bn])
        dofs_changed = np.array(dofs_changed) + f2.num.GLOBAL_dofs \
                       + f3.num.GLOBAL_dofs

        dofs_changed_2 = list() # the rows that been changed in the global matrix in the second apply_strong_BC
        for bn in BNS_com: # all the boundary dofs.
            dofs_changed_2.extend(t2GBD[bn])
        dofs_changed_2 = np.array(dofs_changed_2) + f2.num.GLOBAL_dofs \
                       + f3.num.GLOBAL_dofs

        for i in dofs_changed:
            assert aA[i].nnz == 1 and ab[i] != 0

            if i in dofs_changed_2:
                assert aA[i,i] == 0


        NOT_changed = list() # the rows that not been changed in the global matrix.
        for i in range(aA.shape[0]):
            if i in dofs_changed:
                pass
            else:
                NOT_changed.append(i)
        ab_SUB = ab[NOT_changed]
        assert np.all(ab_SUB==0), f"not change places must be all zero!"
        aA_SUB = aA[NOT_changed,:]
        aA_SUB = aA_SUB[f2.num.GLOBAL_dofs:, f2.num.GLOBAL_dofs:]
        assert aA_SUB.nnz == 0, f"not change places must be all zero!"

        # noinspection PyUnboundLocalVariable
        for i in changed: # the changed rows after the first time `apply_strong_BC`.
            assert aA[i].nnz == 1 and aA[i, i] == 1 and ab[i] != 0


    # check changes in block[2][0] .............................
    GM_F2 = f2.numbering.gathering
    gn_T2 = dict()
    gn_F2 = dict()
    for e in t2BC2:
        gn_T2[e] = GM_t2[e][t2BC2.interpreted_as.local_dofs[e]]
        gn_F2[e] = GM_F2[e][f2BC2.dofs.interpreted_as.local_dofs[e]]
    gn_T2 = cOmm.gather(gn_T2, root=mAster_rank)
    gn_F2 = cOmm.gather(gn_F2, root=mAster_rank)
    if rAnk == mAster_rank:
        AT2 = dict()
        for _ in gn_T2: AT2.update(_)
        AF2 = dict()
        for _ in gn_F2: AF2.update(_)
        for e in AT2:
            gt = AT2[e]
            gf = AF2[e]
            gt = gt + f2.num.GLOBAL_dofs + f3.num.GLOBAL_dofs
            for i, j in zip(gt, gf):
                assert aA[i].nnz == 1 and aA[i,j] == 1 and ab[i] != 0
                assert i not in NOT_changed, f"dof #{i} must be changed"

    return 1


if __name__ == '__main__':
    # mpiexec -n 5 python tests\unittests\tools_.py

    test_TOOLS_NO15_linear_system_apply_BC()


