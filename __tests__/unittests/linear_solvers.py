
import sys
if './' not in sys.path: sys.path.append('./')

from root.config.main import *
import random
from scipy import sparse as spspa
from tools.linear_algebra.data_structures.global_matrix.main import GlobalVector, GlobalMatrix, LocallyFullVector
from tools.linear_algebra.solvers.parallel.GMRES.main import GMRES
from tools.linear_algebra.solvers.parallel.allocator import ParallelSolverDistributor


def ___generate_A_b_of_Manguoglu_Paper___():
    """See [A domain-decomposing parallel sparse linear system solver] by Murat Manguoglu On [Journal of Computational
    and Applied Mathematics]
    """
    A = np.array([(0.2, 1, -1, 0, 0.01, 0, 0, 0, -0.01),
                  (0.01, 0.3, 0, 0, 0, 0, 0, 0, 0),
                  (-0.1, 0, 0.4, 0, 0.3, 0, 0, 0, 0),
                  (0, 0, 0, 0.3, 0.6, 2, 0, 0, 0),
                  (0, -0.2, 0, 0, 0.4, 0, 0, 0, 1.1),
                  (0, 0, 0, -0.2, 0.1, 0.5, 0, 0, 0),
                  (1.2, 0, 0, 0, 0, 0, 0.4, 0.02, 3.0),
                  (0, 0, 0, 0, 0, 0, 2.0, 0.5, 0),
                  (0, 0, 0, 0, 0, 0, 0, 0.1, 0.6)])
    b = np.array([(1,),
                  (1,),
                  (1,),
                  (1,),
                  (1,),
                  (1,),
                  (1,),
                  (1,),
                  (1,)])
    A = spspa.csc_matrix(A)
    b = spspa.csc_matrix(b)

    if rAnk != mAster_rank:
        Ar = spspa.random(9, 9, random.random()/5, format='csc')
        br = spspa.random(9, 1, random.random()/5, format='csc')
    else:
        Ar = spspa.csc_matrix((9,9))
        br = spspa.csc_matrix((9,1))

    Ar0 = cOmm.gather(Ar, root=mAster_rank)
    br0 = cOmm.gather(br, root=mAster_rank)
    if rAnk == mAster_rank:
        AR0 = np.sum(Ar0)
        BR0 = np.sum(br0)
        Ar = A - AR0
        br = b - BR0

    A = GlobalMatrix(spspa.csr_matrix(Ar))
    b = GlobalVector(spspa.csr_matrix(br))

    M = A.___PRIVATE_gather_M_to_core___()
    V = b.___PRIVATE_gather_V_to_core___()
    if rAnk==mAster_rank:
        np.testing.assert_array_almost_equal(M.toarray(),
                                             np.array([(0.2 , 1, -1, 0, 0.01, 0, 0, 0, -0.01),
                                                       (0.01, 0.3, 0, 0, 0, 0, 0, 0, 0),
                                                       (-0.1, 0, 0.4, 0, 0.3, 0, 0, 0, 0),
                                                       (0   , 0, 0, 0.3, 0.6, 2, 0, 0, 0),
                                                       (0   , -0.2, 0, 0, 0.4, 0, 0, 0, 1.1),
                                                       (0   , 0, 0, -0.2, 0.1, 0.5, 0, 0, 0),
                                                       (1.2 , 0, 0, 0, 0, 0, 0.4, 0.02, 3.0),
                                                       (0   , 0, 0, 0, 0, 0, 2.0, 0.5, 0),
                                                       (0   , 0, 0, 0, 0, 0, 0, 0.1, 0.6)]))
        np.testing.assert_array_almost_equal(V, np.array([1,1,1,1,1,1,1,1,1]))

    if sIze > 1: assert not A.IS.master_dominating, f"designed to be."

    return A, b


def test_LinearSolver_No0_GMRES():
    """"""
    if rAnk == mAster_rank:
        print("}}} [test_LinearSolver_No0_GMRES] ...... ", flush=True)
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
            Ar = np.zeros((3,3)) # An even can be empty in some cores.
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
    A = GlobalMatrix(spspa.csc_matrix(Ar))
    b = GlobalVector(spspa.csc_matrix(br))
    X0 = LocallyFullVector(np.zeros((3,)))

    x0, info, beta, ITER, message = ParallelSolverDistributor("GMRES", routine='1', name='GMRES_test_mpi_v2')(
        A, b, X0, restart=3, preconditioner=('Jacobian', dict()), COD=False, plot_residuals=False)
    x0 = x0.V
    np.testing.assert_array_almost_equal(x0, np.array([-2.1810344827586, 1.8362068965517, -0.5948275862068]))

    x0, info, beta, ITER, message = ParallelSolverDistributor("GMRES", routine='0', name='GMRES_test_mpi_v0')(
        A, b, X0, restart=3, preconditioner=('Jacobian', dict()), COD=False, plot_residuals=False)
    x0 = x0.V
    np.testing.assert_array_almost_equal(x0, np.array([-2.1810344827586, 1.8362068965517, -0.5948275862068]))
    assert ITER == 1

    x0, info, beta, ITER, message = GMRES(routine='0', name='GMRES_test_mpi_v0-1')(
        A, b, X0, restart=3, preconditioner=None, COD=False, plot_residuals=False)
    x0 = x0.V
    np.testing.assert_array_almost_equal(x0, np.array([-2.1810344827586, 1.8362068965517, -0.5948275862068]))
    assert ITER == 1



    M = A.___PRIVATE_gather_M_to_core___()
    if rAnk == mAster_rank:
        np.testing.assert_array_almost_equal(M.toarray(),
                                             np.array([(1, 4, 7),
                                                       (2, 9, 7),
                                                       (5, 8, 3)]))

    A, b = ___generate_A_b_of_Manguoglu_Paper___()
    X0 = LocallyFullVector(np.zeros((9,)))

    x0, info, beta, ITER, message = GMRES(routine='1')(A, b, X0, restart=9, preconditioner=None, COD=False)
    x0 = x0.V
    np.testing.assert_array_almost_equal(x0, np.array([-3.23891085,  3.44129703,  1.7765975 , -2.7063454 , -0.11510028,
                                                        0.94048189,  0.36495389,  0.54018445,  1.57663592]))

    x0, info, beta, ITER, message = GMRES(routine='0')(A, b, X0, restart=9, preconditioner=None, COD=False)
    x0 = x0.V
    np.testing.assert_array_almost_equal(x0, np.array([-3.23891085,  3.44129703,  1.7765975 , -2.7063454 , -0.11510028,
                                                        0.94048189,  0.36495389,  0.54018445,  1.57663592]))

    x0, info, beta, ITER, message = ParallelSolverDistributor("GMRES", routine='0')(A, b, X0, restart=9, preconditioner=('Jacobian', dict()),
                                                                       COD=False)
    x0 = x0.V
    np.testing.assert_array_almost_equal(x0, np.array([-3.23891085,  3.44129703,  1.7765975 , -2.7063454 , -0.11510028,
                                                        0.94048189,  0.36495389,  0.54018445,  1.57663592]))

    x0, info, beta, ITER, message = ParallelSolverDistributor("GMRES", routine='auto')(A, b, X0, restart=9, preconditioner=None, COD=False,
                                                                       )
    x0 = x0.V
    np.testing.assert_array_almost_equal(x0, np.array([-3.23891085,  3.44129703,  1.7765975 , -2.7063454 , -0.11510028,
                                                        0.94048189,  0.36495389,  0.54018445,  1.57663592]))

    x0, info, beta, ITER, message = ParallelSolverDistributor("GMRES", routine='auto')(A, b, X0, restart=9, preconditioner=None, COD=False,
                                                                       loading_factor=0) # make sure we use parallel routine
    x0 = x0.V
    np.testing.assert_array_almost_equal(x0, np.array([-3.23891085,  3.44129703,  1.7765975 , -2.7063454 , -0.11510028,
                                                        0.94048189,  0.36495389,  0.54018445,  1.57663592]))


    return 1

def test_LinearSolver_No1_BiCGSTAB():
    """"""
    if rAnk == mAster_rank:
        print("--- [test_LinearSolver_No1_BiCGSTAB] ...... ", flush=True)
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
            Ar = np.zeros((3,3)) # `A` even can be empty in some cores.
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
    A = GlobalMatrix(spspa.csc_matrix(Ar))
    b = GlobalVector(spspa.csc_matrix(br))
    X0 = LocallyFullVector(np.zeros((3,)))

    x0, info, beta, ITER, message = \
    ParallelSolverDistributor("BiCGSTAB")(A, b, X0, maxiter=10, preconditioner=None, COD=False)
    x0 = x0.V
    np.testing.assert_array_almost_equal(x0, np.array([-2.1810344827586, 1.8362068965517, -0.5948275862068]))

    x0, info, beta, ITER, message = \
    ParallelSolverDistributor("BiCGSTAB")(A, b, X0, maxiter=10, preconditioner=('Jacobian', dict()), COD=False)
    x0 = x0.V
    np.testing.assert_array_almost_equal(x0, np.array([-2.1810344827586, 1.8362068965517, -0.5948275862068]))

    x0, info, beta, ITER, message = \
    ParallelSolverDistributor("BiCGSTAB")(A, b, X0, maxiter=10, preconditioner=None, COD=False)
    x0 = x0.V
    np.testing.assert_array_almost_equal(x0, np.array([-2.1810344827586, 1.8362068965517, -0.5948275862068]))


    M = A.___PRIVATE_gather_M_to_core___()
    if rAnk == mAster_rank:
        np.testing.assert_array_almost_equal(M.toarray(),
                                             np.array([(1, 4, 7),
                                                       (2, 9, 7),
                                                       (5, 8, 3)]))

    return 1

def test_LinearSolver_No2_LooseGMRES():
    """"""
    if rAnk == mAster_rank:
        print("--- [test_LinearSolver_No2_LooseGMRES] ...... ", flush=True)

    A, b = ___generate_A_b_of_Manguoglu_Paper___()
    X0 = LocallyFullVector(np.zeros((9,)))

    x0, info, beta, ITer, message = ParallelSolverDistributor("LGMRES")(A, b, X0, m=6, k=2, atol=1e-9,
                                                            maxiter=100, preconditioner=None, COD=False)
    x0 = x0.V
    np.testing.assert_array_almost_equal(x0, np.array([-3.23891085,  3.44129703,  1.7765975 , -2.7063454 , -0.11510028,
                                                        0.94048189,  0.36495389,  0.54018445,  1.57663592]))

    x0, info, beta, ITER, message = ParallelSolverDistributor("LGMRES")(A, b, X0, m=6, k=2, atol=1e-9,
                                                            maxiter=100, preconditioner=('Jacobian', dict()),
                                                            COD=False)
    x0 = x0.V
    np.testing.assert_array_almost_equal(x0, np.array([-3.23891085,  3.44129703,  1.7765975 , -2.7063454 , -0.11510028,
                                                        0.94048189,  0.36495389,  0.54018445,  1.57663592]))

    x0, info, beta, Iter, message = ParallelSolverDistributor("GMRES")(A, b, X0, restart=8, atol=1e-9,
                                                            maxiter=100, preconditioner=None, COD=False)
    assert Iter > ITer > ITER


    return 1

def test_LinearSolver_No3_direct():
    """"""
    if rAnk == mAster_rank:
        print("ddd [test_LinearSolver_No3_direct] ...... ", flush=True)


    A, b = ___generate_A_b_of_Manguoglu_Paper___()

    x0, info, beta, ITer, message = ParallelSolverDistributor("direct")(A, b, COD=False)
    x0 = x0.V
    np.testing.assert_array_almost_equal(x0, np.array([-3.23891085,  3.44129703,  1.7765975 , -2.7063454 ,
                                                       -0.11510028,
                                                        0.94048189,  0.36495389,  0.54018445,  1.57663592]))

    np.testing.assert_almost_equal(A.condition.condition_number, 85.3100212781)

    x0, info, beta, ITer, message = ParallelSolverDistributor("direct")(A, b, COD=True)
    x0 = x0.V
    np.testing.assert_array_almost_equal(x0, np.array([-3.23891085,  3.44129703,  1.7765975 , -2.7063454 ,
                                                       -0.11510028,
                                                        0.94048189,  0.36495389,  0.54018445,  1.57663592]))

    assert A.IS.master_dominating # the COD=True has triggered this!
    np.testing.assert_almost_equal(A.condition.condition_number, 85.3100212781)

    return 1











if __name__ == '__main__':
    # mpiexec -n 4 python tests\unittests\linear_solvers.py

    test_LinearSolver_No2_LooseGMRES()
    test_LinearSolver_No1_BiCGSTAB()
    test_LinearSolver_No0_GMRES()
    test_LinearSolver_No3_direct()