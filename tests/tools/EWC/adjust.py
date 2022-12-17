# -*- coding: utf-8 -*-
import sys
if './' not in sys.path:
    sys.path.append('./')

from tools.elementwiseCache.dataStructures.objects.sparseMatrix.main import EWC_SparseMatrix
from tools.elementwiseCache.dataStructures.objects.columnVector.main import EWC_ColumnVector
from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller
from root.config.main import RANK, MASTER_RANK, COMM, np
import random
from scipy.sparse import random as sci_rand


def test_EWC_adjust():
    """"""
    if RANK == MASTER_RANK:
        print("ADJ [test_EWC_adjust] ...... ", flush=True)
        i, j, k = random.randint(1, 3), random.randint(1, 3), random.randint(1, 3)
        a, b, c = random.randint(1, 3), random.randint(1, 3), random.randint(1, 3)
        Dirichlet_boundaries = random.sample(
            ['North', 'West', 'South', 'Back', 'Front', "East"], random.randint(0, 6))
    else:
        i, j, k, a, b, c = [None for _ in range(6)]
        Dirichlet_boundaries = None

    i, j, k, a, b, c = COMM.bcast([i, j, k, a, b, c], root=MASTER_RANK)
    Dirichlet_boundaries = COMM.bcast(Dirichlet_boundaries, root=MASTER_RANK)

    mesh = MeshGenerator('crazy', c=0.1)([i, j, k])
    Neumann_boundaries = list()
    BNS = mesh.boundaries.names
    for bn in BNS:
        if bn not in Dirichlet_boundaries:
            Neumann_boundaries.append(bn)

    space = SpaceInvoker('polynomials')([a, b, c])
    FC = FormCaller(mesh, space)

    def u(t, x, y, z): return 0 * x * y * z * t
    def v(t, x, y, z): return 0 * x * y * z * t
    def w(t, x, y, z): return 0 * x * y * z * t

    velocity = FC('vector', (u, v, w))

    f1 = FC('1-f', hybrid=True)
    t1 = FC('1-adt')

    f1.BC.CF = velocity
    velocity.current_time = 0
    f1.BC.boundaries = Neumann_boundaries
    t1.BC.boundaries = Neumann_boundaries
    col_pc = f1.BC.interpret
    row_pd = t1.BC.interpret

    T = dict()
    for i in mesh.elements:
        T[i] = sci_rand(t1.num.basis, f1.num.basis, density=1, format='csr')
    T = EWC_SparseMatrix(mesh, T)
    T.gathering_matrices = (t1, f1)
    T = T.adjust.identify_rows_according_to(row_pd, col_pc)

    b = dict()
    for i in mesh.elements:
        b[i] = sci_rand(t1.num.basis, 1, density=1, format='csc')
    b = EWC_ColumnVector(mesh, b)
    b.gathering_matrix = t1
    b = b.adjust.set_entries_according_to(row_pd, col_pc)

    Ta = T.assembled
    ba = b.assembled
    TT = Ta.do.gather_M_to_core()
    bb = ba.do.gather_V_to_core()

    if RANK == MASTER_RANK:
        for _, Ti in enumerate(TT):
            if Ti.nnz == 1:
                assert Ti[0, Ti.indices[0]] == 1
                assert bb[_] == 0

    T = t1.matrices.trace
    NT = T.adjust.identify_rows_according_to(row_pd, col_pc)
    for i in T:
        assert np.all(T[i].indptr == NT[i].indptr)

    T = T.adjust.identify_rows_according_to(row_pd, col_pc)
    D = EWC_SparseMatrix(mesh, (t1.num.basis, t1.num.basis))
    D.gathering_matrices = (t1, t1)
    T.gathering_matrices = (t1, f1)

    t1.prime.BC.CF = velocity.components.T_perp
    t1.prime.BC.CF.current_time = 0
    t1.BC.boundaries = Dirichlet_boundaries
    tpc = t1.BC.interpret
    D = D.adjust.identify_rows_according_to(tpc)
    T = T.adjust.clear_rows_according_to(tpc)
    b = b.adjust.set_entries_according_to(tpc, tpc)

    Da = D.assembled
    Ta = T.assembled
    ba = b.assembled
    DD = Da.do.gather_M_to_core()
    TT = Ta.do.gather_M_to_core()
    bb = ba.do.gather_V_to_core()

    if RANK == MASTER_RANK:
        for _, Ti in enumerate(TT):
            if Ti.nnz == 1:
                assert Ti[0, Ti.indices[0]] == 1
                assert bb[_] == 0
            elif Ti.nnz == 2:
                assert np.sum(Ti) == 0
            else:
                assert Ti.nnz == 0
                assert bb[_] == 0
                Di = DD[_]
                assert Di.nnz == 1
                assert Di[0, _] == 1

    return 1


if __name__ == '__main__':
    # mpiexec -n 4 python __tests__/unittests/EWC/adjust.py
    test_EWC_adjust()
