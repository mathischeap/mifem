# -*- coding: utf-8 -*-
import sys
if './' not in sys.path: sys.path.append('./')


from root.config.main import *
import random
from objects.CSCG._3d.master import SpaceInvoker, FormCaller
from tools.elementwiseCache.dataStructures.objects.sparseMatrix.main import EWC_ColumnVector
from tools.elementwiseCache.dataStructures.operators.concatenate.main import concatenate
from scipy import sparse as spspa
from tests.objects.CSCG._3d.randObj.mesh import random_mesh_of_elements_around as _3d_RANDOM_MESH_


def test_LinearAlgebra_EWC_No0_ColumnVector():
    """"""
    if RANK == MASTER_RANK:
        print("&&& [test_LinearAlgebra_No0_EWC_ColumnVector] ...... ", flush=True)

    if RANK == MASTER_RANK:
        load = random.randint(10,100)
        IH = [True, False][random.randint(0,1)]
        l = random.randint(2, 3)
        m = random.randint(1, 3)
        n = random.randint(2, 4)
    else:
        load = None
        IH = None
        l, m, n = None, None, None
    load, IH = COMM.bcast([load, IH], root=MASTER_RANK)
    l, m, n = COMM.bcast([l, m, n], root=MASTER_RANK)
    mesh = _3d_RANDOM_MESH_(load)
    space = SpaceInvoker('polynomials')([('Lobatto', l), ('Lobatto', m), ('Lobatto', n)])
    FC = FormCaller(mesh, space)

    def u(t, x, y, z): return np.cos(np.pi*x) + np.sin(np.pi*y) * np.sin(np.pi*z-0.125)**2 + t/2
    def v(t, x, y, z): return np.sin(np.pi*x) + np.sin(np.pi*y) * np.sin(np.pi*z-0.125)**2 + t/2
    def w(t, x, y, z): return np.sin(np.pi*x) + np.cos(np.pi*y) * np.cos(np.pi*z-0.125)**2 + t
    def p(t, x, y, z): return np.cos(np.pi*x) + np.sin(np.pi*y) * np.sin(np.pi*z-0.125)**2 + t/2

    scalar = FC('scalar', p)
    vector = FC('vector', (u,v,w))
    f0 = FC('0-f', is_hybrid=False)
    f1 = FC('1-f', is_hybrid=False)
    f0.CF = scalar
    scalar.current_time = 0
    f0.discretize()
    f1.CF = vector
    vector.current_time = 0
    f1.discretize()

    v0 = f0.cochain.EWC
    v1 = f1.cochain.EWC
    v2 = EWC_ColumnVector(mesh, 100)
    V = concatenate([v0, v1, v2])

    for i in V:
        v0_ = v0[i]
        v1_ = v1[i]
        v2_ = v2[i]
        V_ = V[i]
        vvv = spspa.vstack([v0_, v1_, v2_], format='csc')
        ZERO = V_ - vvv
        assert ZERO.nnz == 0


    #-------- tests of customize ------------------------------------------------------------------
    V = concatenate([v0, v1])
    GLOBAL_num_dofs = V.gathering_matrix.global_num_dofs
    if RANK == MASTER_RANK:
        num_samples = random.randint(1,5)
        indices = random.sample(range(0, GLOBAL_num_dofs), num_samples)
        values = [random.random() for __ in range(num_samples)]
        renew_values = [random.random() for __ in range(num_samples)]

    else:
        indices = None
        values = None
        renew_values = None
    indices, values, renew_values = COMM.bcast([indices, values, renew_values], root=MASTER_RANK)
    for i, v in zip(indices, values):
        V.customize.set_assembled_V_i_to(i, v)

    Va = V.assembled
    Va = Va.do.gather_V_to_core()
    if RANK == MASTER_RANK:
        np.testing.assert_array_almost_equal(Va[indices] - values, 0)

    for i, v in zip(indices, renew_values):
        V.customize.set_assembled_V_i_to(i, v)

    Va = V.assembled
    Va = Va.do.gather_V_to_core()
    if RANK == MASTER_RANK:
        np.testing.assert_array_almost_equal(Va[indices] - renew_values, 0)


    return 1





if __name__ == '__main__':
    # mpiexec -n 4 python __tests__/unittests/EWC/column_vector.py

    test_LinearAlgebra_EWC_No0_ColumnVector()
