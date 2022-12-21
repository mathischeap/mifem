# -*- coding: utf-8 -*-
"""
@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
import sys
if './' not in sys.path:
    sys.path.append('./')

from root.config.main import *
import os
from time import sleep
from scipy import sparse as spspa
from tools.iterators.simple import SimpleIterator
import random
from tools.miLinearAlgebra.dataStructures.globalMatrix.main import GlobalMatrix
from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller
from objects.CSCG._2d.master import MeshGenerator as MeshGenerator2D
from objects.CSCG._2d.master import SpaceInvoker as SpaceInvoker2D
from objects.CSCG._2d.master import FormCaller as FormCaller2D
from tools.elementwiseCache.gathering.regular.chain_matrix.main import Chain_Gathering_Matrix
from tools.elementwiseCache.dataStructures.operators.concatenate.main import concatenate
from tools.elementwiseCache.dataStructures.operators.bmat.main import bmat
from tools.miLinearAlgebra.linearSystem.main import LinearSystem
from tools.elementwiseCache.dataStructures.objects.sparseMatrix.main import EWC_ColumnVector
from tools.run.reader import ParallelMatrix3dInputRunner, RunnerDataReader
from tests.objects.CSCG._3d.randObj.form_caller import random_mesh_and_space_of_total_load_around


def ___TEST_SOLVER___(tk, tk1):
    """
    Parameters
    ----------
    tk :
    tk1 ：

    Returns
    -------
    exit_code : The standard exit code.
    shut_down : If it is ``True``, the outer iterator will `shutdown` immediately.
    message : The solver message.
    output1 :
    output2 :
    """
    sleep(0.01)
    assert tk1 > tk
    if RANK == MASTER_RANK:
        ot1 = (tk+tk1)/2
        ot2 = (tk+0.1*tk1)/3
    else:
        ot1 = None
        ot2 = None
    ot1, ot2 = COMM.bcast([ot1, ot2], root=MASTER_RANK)
    return 1, 0, 'random message', ot1, ot2


def test_TOOLS_NO1_iterator():
    if RANK == MASTER_RANK:
        print("=== [test_TOOLS_NO1_iterator] ...... ", flush=True)
    SI = SimpleIterator(
        t0=0, dt=0.1, max_steps=10,
        auto_save_frequency=5,
        monitor_factor=0,
        RDF_filename='RDF_filename',
        save_to_mitr=True,
        name=None
    )
    SI(___TEST_SOLVER___, [0, 0])
    SI.run()

    S2 = SimpleIterator.read("RDF_filename.mitr")

    if RANK == MASTER_RANK:
        assert S2._solver_dir_ == SI._solver_dir_
        assert S2._solver_source_code_ == SI._solver_source_code_
        np.testing.assert_array_equal(SI.RDF.to_numpy(), S2.RDF.to_numpy())

        os.remove('RDF_filename.csv')
        os.remove('RDF_filename.mitr')
    else:
        assert S2 is None, "Iterator only read to master core. It is only for retrieve info, not for resume simulation."
    return 1


def test_TOOLS_NO6_send_GM_in_parts_test():
    if RANK == MASTER_RANK:
        print("=== [test_TOOLS_NO6_send_GM_in_parts_test] ....", flush=True)

    for T in range(50):
        if RANK == MASTER_RANK:
            i = random.randint(0, 1)
            j = random.randint(10, 99)
            k = random.randint(100, 199)
        else:
            i, j, k = None, None, None
        i = COMM.bcast(i, root=MASTER_RANK)
        j = COMM.bcast(j, root=MASTER_RANK)
        k = COMM.bcast(k, root=MASTER_RANK)
        L = random.randint(2, 4)
        empty = sorted(random.sample(range(0, j), int(j / L)))

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
            A = A.tocsc()
        else:
            raise Exception()
        A0 = COMM.gather(A, root=MASTER_RANK)
        if RANK == MASTER_RANK:
            A0 = np.sum(A0)

        A = GlobalMatrix(A)

        A = A.do.gather_M_to_core(clean_local=True, splitting_factor=k)
        if RANK == MASTER_RANK:
            np.testing.assert_array_almost_equal(A.toarray(), A0.toarray())

    return 1


def test_TOOLS_NO7_linear_algebra_EWC_test():
    if RANK == MASTER_RANK:
        print("::: [test_TOOLS_NO7_linear_algebra_EWC_test] ....", flush=True)

    for _ in range(2):
        if RANK == MASTER_RANK:
            c = random.uniform(0, 0.15)
            i = random.randint(2, 3)
            j = random.randint(1, 3)
            k = random.randint(2, 4)
            L = random.randint(2, 3)
            m = random.randint(1, 3)
            n = random.randint(2, 4)
        else:
            c, i, j, k, L, m, n = None, None, None, None, None, None, None
        c, i, j, k, L, m, n = COMM.bcast([c, i, j, k, L, m, n], root=MASTER_RANK)
        if c < 0.1:
            c = 0
        mesh = MeshGenerator('crazy', c=c)([f'Lobatto:{i}', f'Lobatto:{j}', f'Lobatto:{k}'], EDM='debug')
        space = SpaceInvoker('polynomials')([('Lobatto', L), ('Lobatto', m), ('Lobatto', n)])
        FC = FormCaller(mesh, space)

        f0 = FC('0-f', hybrid=False)
        f1 = FC('1-f', hybrid=False)
        f2 = FC('2-f', hybrid=False)

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

        blocks = ([M2, E21],
                  [None, M1])

        BMAT = bmat(blocks)

        for i in BMAT:
            bi = BMAT[i]
            m2 = M2[i]
            e21 = E21[i]
            m1 = M1[i]
            Bi = spspa.bmat(([m2, e21], [None, m1]), format='csc')
            np.testing.assert_almost_equal(np.sum(np.abs(bi - Bi)), 0)

    return 1


def test_TOOLS_NO8_GlobalMatrix_dot_product_test():
    """"""
    if RANK == MASTER_RANK:
        print("... [test_TOOLS_NO8_GlobalMatrix_dot_product_test] ....", flush=True)

    for _ in range(20):  # do multiple random tests.
        # A @ B @ C, A: (i,j), B: (j, k), C:(k, l)
        if RANK == MASTER_RANK:
            i = random.randint(1, SIZE * 2)
            j = random.randint(1, SIZE * 2)
            k = random.randint(1, SIZE * 2)
            L = random.randint(1, SIZE * 2)

            type_A = random.randint(0, 1)
            type_B = random.randint(0, 1)
            type_C = random.randint(0, 1)
        else:
            i, j, k, L = None, None, None, None
            type_A, type_B, type_C = None, None, None
        i, j, k, L = COMM.bcast([i, j, k, L], root=MASTER_RANK)
        type_A, type_B, type_C = COMM.bcast([type_A, type_B, type_C], root=MASTER_RANK)
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
            C = ___generate_random_row_major_GM___(k, L)
        else:
            C = ___generate_random_col_major_GM___(k, L)

        ABC = A @ B @ C
        a = A.do.gather_M_to_core()
        b = B.do.gather_M_to_core()
        c = C.do.gather_M_to_core()
        abc = ABC.do.gather_M_to_core()

        if RANK == MASTER_RANK:
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
        if s < 0.02:
            s = 0
    if RANK == MASTER_RANK:

        random_list = random.sample(range(0, i), i)
        distribution = [i // SIZE + (1 if x < i % SIZE else 0) for x in range(SIZE)]
        no_empty_rows = list()
        _ = 0
        for r in range(SIZE):
            no_empty_rows.append(random_list[_:_+distribution[r]])
            _ += distribution[r]

    else:
        no_empty_rows = None

    no_empty_rows = COMM.scatter(no_empty_rows, root=MASTER_RANK)

    _ = spspa.random(i, j, s, format='csr')

    A = spspa.lil_matrix((i, j))
    A[no_empty_rows, :] = _[no_empty_rows, :]
    A = A.tocsr()
    A = GlobalMatrix(A)
    A.whether.regularly_distributed = 'row'
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
        if s < 0.02:
            s = 0
    if RANK == MASTER_RANK:

        random_list = random.sample(range(0, j), j)
        distribution = [j // SIZE + (1 if x < j % SIZE else 0) for x in range(SIZE)]
        no_empty_cols = list()
        _ = 0
        for r in range(SIZE):
            no_empty_cols.append(random_list[_:_+distribution[r]])
            _ += distribution[r]

    else:
        no_empty_cols = None

    no_empty_cols = COMM.scatter(no_empty_cols, root=MASTER_RANK)

    _ = spspa.random(i, j, s, format='csc')

    A = spspa.lil_matrix((i, j))
    A[:, no_empty_cols] = _[:, no_empty_cols]
    A = A.tocsc()
    A = GlobalMatrix(A)
    A.whether.regularly_distributed = 'column'
    A.___PRIVATE_self_regularity_checker___()

    return A


def test_TOOLS_NO9_test_Chained_Gathering_Matrix():
    """"""
    if RANK == MASTER_RANK:
        print("+++ [test_TOOLS_NO9_test_Chained_Gathering_Matrix] ....", flush=True)

    for _ in range(3):
        if RANK == MASTER_RANK:
            i = random.randint(2, 3)
            j = random.randint(1, 3)
            k = random.randint(1, 2)
            L = random.randint(2, 3)
            m = random.randint(1, 2)
            n = random.randint(2, 3)

            mid = random.randint(0, 1)
        else:
            i, j, k, L, m, n = None, None, None, None, None, None
            mid = None
        i, j, k, L, m, n, mid = COMM.bcast([i, j, k, L, m, n, mid], root=MASTER_RANK)

        MID = ('crazy', 'bridge_arch_cracked')[mid]

        mesh = MeshGenerator(MID)([f'Lobatto:{i}', f'Lobatto:{j}', f'Lobatto:{k}'], EDM='debug')
        space = SpaceInvoker('polynomials')([('Lobatto', L), ('Lobatto', m), ('Lobatto', n)])
        FC = FormCaller(mesh, space)
        f0 = FC('0-f', hybrid=False)
        f1 = FC('1-f', hybrid=False)
        f2 = FC('2-f', hybrid=False)
        f3 = FC('3-f', hybrid=False)

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

        NUM0 = GM0.global_num_dofs
        NUM1 = GM1.global_num_dofs
        for i in CGM1:
            cgm = CGM1[i]
            gm0 = GM0[i].full_vector
            gm1 = GM1[i].full_vector
            gm2 = GM2[i].full_vector
            c_g_m = np.concatenate([gm0, gm1+NUM0, gm2+NUM0+NUM1])
            np.testing.assert_array_equal(cgm, c_g_m)

        assert CGM2.global_num_dofs == GM1.global_num_dofs + GM2.global_num_dofs + GM3.global_num_dofs

        CGM3 = Chain_Gathering_Matrix([GM0, GM1, GM2, GM3])
        CGM4 = Chain_Gathering_Matrix([GM0, GM0, GM2, GM3])
        CGM5 = Chain_Gathering_Matrix([GM0, GM3])
        CGM6 = Chain_Gathering_Matrix([GM2, ])

        assert CGM0.mesh_type == '_3dCSCG'
        assert CGM1.mesh_type == '_3dCSCG'
        assert CGM2.mesh_type == '_3dCSCG'
        assert CGM3.mesh_type == '_3dCSCG'
        assert CGM4.mesh_type == '_3dCSCG'
        assert CGM5.mesh_type == '_3dCSCG'
        assert CGM6.mesh_type == '_3dCSCG'

        assert CGM6.global_num_dofs == GM2.global_num_dofs
        assert CGM0.global_num_dofs == 2 * GM1.global_num_dofs
        assert CGM4.global_num_dofs == 2 * GM0.global_num_dofs + GM2.global_num_dofs + GM3.global_num_dofs

        if RANK == MASTER_RANK:
            COMM_CHECK = random.randint(3, 5)
            i = random.randint(0, 6)
        else:
            COMM_CHECK = None
            i = None
        i, COMM_CHECK = COMM.bcast([i, COMM_CHECK], root=MASTER_RANK)
        CGM = (CGM0, CGM1, CGM2, CGM3, CGM4, CGM5, CGM6)[i]
        assert CGM.mesh_type == '_3dCSCG'

        gnd_S = list()
        for gm in CGM.GMs:
            gnd_S.append(gm.global_num_dofs)

        for i in CGM:
            cgm = CGM[i]
            assert len(set(cgm)) == len(cgm), "Having repeated numbering? BAD BAD BAD!"
            gm_list = list()
            for j, gm in enumerate(CGM.GMs):
                gm_list.append(gm[i].full_vector+sum(gnd_S[0:j]))
            c_g_m = np.concatenate(gm_list)
            np.testing.assert_array_equal(cgm, c_g_m)

        GND = CGM.global_num_dofs

        n = 0

        for m in range(GND):
            _I = CGM.do.find.elements_and_local_indices_of_dof(m)
            if _I is None:
                exclude_list = list()
            else:
                _I = _I[0]
                assert isinstance(_I, list)
                exclude_list = _I
                for i in _I:
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
                GET = COMM.gather(GET_it, root=MASTER_RANK)
                if RANK == MASTER_RANK:
                    assert sum(GET) >= 1, f"At least in one core we will find element(s) contain dof {m}!"

    return 1


def test_TOOLS_NO10_test_EWC_SparseMatrix_Customize():
    """"""
    if RANK == MASTER_RANK:
        print("||| [test_TOOLS_NO10_test_EWC_SparseMatrix_Customize] ....", flush=True)

    for _ in range(1):
        if RANK == MASTER_RANK:
            i = random.randint(2, 3)
            j = random.randint(1, 3)
            k = random.randint(1, 2)
            L = random.randint(2, 3)
            m = random.randint(1, 2)
            n = random.randint(2, 3)

            mid = random.randint(0, 1)
            ish = random.randint(0, 1)
        else:
            i, j, k, L, m, n = None, None, None, None, None, None
            mid = None
            ish = None
        i, j, k, L, m, n, mid, ish = COMM.bcast([i, j, k, L, m, n, mid, ish], root=MASTER_RANK)

        MID = ('crazy', 'bridge_arch_cracked')[mid]
        ISH = (True, False)[ish]

        mesh = MeshGenerator(MID)([f'Lobatto:{i}', f'Lobatto:{j}', f'Lobatto:{k}'], EDM='debug')
        space = SpaceInvoker('polynomials')([('Lobatto', L), ('Lobatto', m), ('Lobatto', n)])
        FC = FormCaller(mesh, space)

        f2 = FC('2-f', hybrid=ISH)
        f3 = FC('3-f', hybrid=ISH)

        M3 = f3.matrices.mass
        E32 = f2.matrices.incidence

        f2_GND = f2.num.global_dofs
        f3_GND = f3.num.global_dofs
        GND = f2_GND + f3_GND
        BMAT = [[M3, E32], [E32.T, None]]
        SYSTEM = bmat(BMAT)
        SYSTEM.gathering_matrices = ([f3, f2], [f3, f2])

        if RANK == MASTER_RANK:
            AAA = int(GND / 30)
            if AAA > 10:
                HOW_MANY_LINES = 10
            else:
                HOW_MANY_LINES = AAA
            III = random.sample(range(0, GND), HOW_MANY_LINES)

        else:
            III = None
        III = COMM.bcast(III, root=MASTER_RANK)
        for iii in III:
            SYSTEM.customize.clear_global_row(iii)
        SYSTEM_ASSEMBLED = SYSTEM.assembled

        for jjj in III:
            assert SYSTEM_ASSEMBLED.M[jjj].nnz == 0
        M = SYSTEM_ASSEMBLED.do.gather_M_to_core()
        for jjj in range(GND):
            if jjj not in III:
                if RANK == MASTER_RANK:
                    assert M[jjj].nnz != 0

        SYSTEM = M3
        GND = f3_GND
        if RANK == MASTER_RANK:
            AAA = int(GND / 20)
            if AAA > 10:
                HOW_MANY_LINES = 10
            else:
                HOW_MANY_LINES = AAA
            III = random.sample(range(0, GND), HOW_MANY_LINES)
        else:
            III = None

        III = COMM.bcast(III, root=MASTER_RANK)
        for iii in III:
            SYSTEM.customize.clear_global_row(iii)
        SYSTEM_ASSEMBLED = SYSTEM.assembled

        for jjj in III:
            assert SYSTEM_ASSEMBLED.M[jjj].nnz == 0
        M = SYSTEM_ASSEMBLED.do.gather_M_to_core()
        for jjj in range(GND):
            if jjj not in III:
                if RANK == MASTER_RANK:
                    assert M[jjj].nnz != 0

        f0 = FC('0-f', hybrid=False)
        f1 = FC('1-f', hybrid=False)

        M1 = f1.matrices.mass
        E10 = f0.matrices.incidence

        f0_GND = f0.num.global_dofs
        f1_GND = f1.num.global_dofs
        GND = f0_GND + f1_GND
        BMAT = [[M1, E10], [E10.T, None]]
        SYSTEM = bmat(BMAT)
        SYSTEM.gathering_matrices = ([f1, f0], [f1, f0])

        if RANK == MASTER_RANK:
            AAA = int(GND / 30)
            if AAA > 10:
                HOW_MANY_LINES = 10
            else:
                HOW_MANY_LINES = AAA
            III = random.sample(range(0, GND), HOW_MANY_LINES)
        else:
            III = None
        III = COMM.bcast(III, root=MASTER_RANK)
        for iii in III:
            SYSTEM.customize.clear_global_row(iii)
        SYSTEM_ASSEMBLED = SYSTEM.assembled

        for jjj in III:
            assert SYSTEM_ASSEMBLED.M[jjj].nnz == 0
        M = SYSTEM_ASSEMBLED.do.gather_M_to_core()
        for jjj in range(GND):
            if jjj not in III:
                if RANK == MASTER_RANK:
                    assert M[jjj].nnz != 0

        SYSTEM = bmat(BMAT)
        SYSTEM.gathering_matrices = ([f1, f0], [f1, f0])
        if RANK == MASTER_RANK:
            AAA = int(GND / 30)
            if AAA > 10:
                HOW_MANY_LINES = 10
            else:
                HOW_MANY_LINES = AAA
            III = random.sample(range(0, GND), HOW_MANY_LINES)
            VVV = np.random.rand(HOW_MANY_LINES)
        else:
            III, VVV = None, None
        III, VVV = COMM.bcast([III, VVV], root=MASTER_RANK)  # we will set the value at M[i,j], i in III, j in JJJ.

        for i, v in zip(III, VVV):
            SYSTEM.customize.set_assembled_M_ij_to(i, i, v)
        SYSTEM_ASSEMBLED = SYSTEM.assembled
        M = SYSTEM_ASSEMBLED.do.gather_M_to_core()

        if RANK == MASTER_RANK:
            for i, v in zip(III, VVV):
                assert M[i, i] == v

        SYSTEM = bmat(BMAT)
        SYSTEM.gathering_matrices = ([f1, f0], [f1, f0])
        if RANK == MASTER_RANK:
            AAA = int(GND / 30)
            if AAA > 10:
                HOW_MANY_LINES = 10
            else:
                HOW_MANY_LINES = AAA
            III = random.sample(range(0, GND), HOW_MANY_LINES)
        else:
            III, VVV = None, None
        III = COMM.bcast(III, root=MASTER_RANK)  # we will set the value at M[i,j], i in III, j in JJJ.

        for i in III:
            SYSTEM.customize.identify_global_row(i)
        SYSTEM_ASSEMBLED = SYSTEM.assembled
        M = SYSTEM_ASSEMBLED.do.gather_M_to_core()

        if RANK == MASTER_RANK:
            for i in III:
                assert M[i].nnz == 1
                assert M[i, i] == 1

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
    if RANK == MASTER_RANK:
        result = np.array(SI.RDF.to_numpy()[-1, :][-2:], dtype=float)
        return result[0]+A, result[1]-A


def test_TOOLS_NO11_test_ParallelMatrix3dInputRunner():
    """"""
    if RANK == MASTER_RANK:
        print("->- [test_TOOLS_NO11_test_ParallelMatrix3dInputRunner] ....", flush=True)

        i1 = random.randint(2, 3)
        i2 = random.randint(1, 4)
        J = random.randint(1, 3)

        T0 = (random.sample(range(0, 2*i1), i1),
              random.sample(range(2*i1+1, 2*i1+1+2*i2), i2))
        DT = ([(random.random()+0.1) * np.pi / 10 for _ in range(i1)],
              [(random.random()+0.1) * np.sqrt(2) / 10 for _ in range(i2)])
        STEPS = random.sample(range(1, 2*J), J)

        A = random.uniform(0, 10)
    else:
        T0, DT, STEPS, A = None, None, None, None

    T0, DT, STEPS, A = COMM.bcast([T0, DT, STEPS, A], root=MASTER_RANK)

    PR1 = ParallelMatrix3dInputRunner(___runner_test_function_11___)
    PR1.iterate(T0, DT, STEPS, writeto='pmr_test.txt', A=A)

    if RANK == MASTER_RANK:

        DR1 = len(PR1._SR_.rdf)
        if DR1 < 5:   # if we have less than 5 runs, we re-run all.
            DEL_ROWS = DR1
        elif DR1 > 15:
            DEL_ROWS = 0   # re-run none.
        elif DR1 > 11:
            DEL_ROWS = random.randint(0, int(DR1/3))
        else:
            DEL_ROWS = random.randint(0, DR1)

        with open('pmr_test.txt', 'r') as f:
            contents = f.readlines()
            LEN1 = len(contents)
            ROWS = random.sample(range(0, LEN1-1), random.randint(1, 10))
            ROWS.extend([-1, -2, -3])
            TO_BE_CHECK = list()
            for r in ROWS:
                TO_BE_CHECK.append(contents[r])

        if DEL_ROWS != 0:
            with open('pmr_test.txt', 'w') as f:
                for con in contents[:-DEL_ROWS]:
                    f.write(con)

    PR2 = ParallelMatrix3dInputRunner(___runner_test_function_11___)
    PR2.iterate(T0, DT, STEPS, writeto='pmr_test.txt', A=A)

    if RANK == MASTER_RANK:
        with open('pmr_test.txt', 'r') as f:
            contents = f.readlines()
            for i, r in enumerate(ROWS):

                if len(TO_BE_CHECK[i]) >= 138:
                    assert TO_BE_CHECK[i][:138] == contents[r][:138]
                else:
                    assert TO_BE_CHECK[i] == contents[r]

    PR3 = RunnerDataReader('pmr_test.txt')

    if RANK == MASTER_RANK:
        D = PR2._SR_.rdf.to_numpy()
    else:
        D = None

    D = COMM.bcast(D, root=MASTER_RANK)
    assert np.all(D == PR3.results.to_numpy())
    assert PR3.___lock_iterate___

    if RANK == MASTER_RANK:
        os.remove('pmr_test.txt')

    return 1


def test_TOOLS_NO12_EWC_assembling_test():
    """"""
    if RANK == MASTER_RANK:
        print("AAA [test_TOOLS_NO12_EWC_assembling_test] ....", flush=True)

    if RANK == MASTER_RANK:
        el1 = random.randint(1, 3)
        el2 = random.randint(1, 3)
        el3 = random.randint(1,  3)
        c = random.uniform(0.0, 0.3)
        if c < 0.1:
            c = 0
    else:
        el1, el2, el3, c = [None for _ in range(4)]
    el1, el2, el3, c = COMM.bcast([el1, el2, el3, c], root=MASTER_RANK)

    mesh = MeshGenerator('crazy', c=c)([el1, el2, el3])
    space = SpaceInvoker('polynomials')([('Lobatto', el3), ('Lobatto', el2), ('Lobatto', el1)])
    FC = FormCaller(mesh, space)
    f0 = FC('0-f', hybrid=False)
    f1 = FC('1-f', hybrid=False)
    f2 = FC('2-f', hybrid=False)
    f3 = FC('2-f', hybrid=False)

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

    f0.CF = scalar
    scalar.current_time = 0
    f0.discretize()
    c0 = f0.cochain.local
    aC0 = ___brutal_force_EWC_vector_assembling___(c0, GM0)
    C0 = f0.cochain.EWC
    AC0 = C0.assembled
    np.testing.assert_array_almost_equal(aC0.toarray(), AC0.V.toarray())

    mesh = MeshGenerator2D('crazy', c=c)([el1, el2])
    space = SpaceInvoker2D('polynomials')([('Lobatto', el3), ('Lobatto', el2)])
    FC = FormCaller2D(mesh, space)
    f0 = FC('0-f-i', hybrid=False)
    f1 = FC('1-f-i', hybrid=False)
    f2 = FC('2-f-i', hybrid=False)
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
    M = bmat(([M0, ], [E10, ]))

    aM = ___brutal_force_EWC_matrix_assembling___(M, *M.gathering_matrices)
    AM = M.assembled
    np.testing.assert_array_almost_equal(aM.toarray(), AM.M.toarray())

    def p(t, x, y): return -0.001 + x + 2.1254*y + t
    scalar = FC('scalar', p)

    f0.CF = scalar
    f0.CF.current_time = 0
    f0.discretize()
    c0 = f0.cochain.local
    aC0 = ___brutal_force_EWC_vector_assembling___(c0, GM0)
    C0 = f0.cochain.EWC
    AC0 = C0.assembled
    np.testing.assert_array_almost_equal(aC0.toarray(), AC0.V.toarray())

    return 1


def ___brutal_force_EWC_matrix_assembling___(EWC, GM0, GM1):
    """Do the brutal force assembling."""
    N0 = GM0.global_num_dofs
    N1 = GM1.global_num_dofs

    AM = spspa.lil_matrix((N0, N1))
    for i in EWC:
        Mi = EWC[i]
        gm0 = GM0[i]
        gm1 = GM1[i]

        for j, m in enumerate(gm0):
            for k, n in enumerate(gm1):
                AM[m, n] += Mi[j, k]

    return AM


def ___brutal_force_EWC_vector_assembling___(EWC, GM):
    """Do the brutal force assembling."""
    N = GM.global_num_dofs

    A = spspa.lil_matrix((N, 1))
    for i in EWC:
        Vi = EWC[i]
        gm = GM[i]

        for j, m in enumerate(gm):
            A[m, 0] += Vi[j]

    return A


def test_TOOLS_NO13_EWC_Customize_CSCG_partial_dofs():
    """"""
    if RANK == MASTER_RANK:
        print("-p- [test_TOOLS_NO13_EWC_Customize_CSCG_partial_dofs] ....", flush=True)

    if RANK == MASTER_RANK:
        ISH = (True, False, False)[random.randint(0, 2)]
        LOAD = random.randint(25, 75)
    else:
        ISH = None
        LOAD = None
    ISH, LOAD = COMM.bcast([ISH, LOAD], root=MASTER_RANK)

    mesh, space = random_mesh_and_space_of_total_load_around(
        LOAD,
        exclude_periodic=True,
        domain_boundary_distribution_regularities='Regular:interfaces-not-shared-by-regions'
    )

    # ------- test partial_dofs only having boundary dofs involved====
    bns = mesh.boundaries.names
    if RANK == MASTER_RANK:
        BNS = list()
        for _ in range(4):
            i = random.randint(1, len(bns))  # how many boundaries to be included.
            BNS.append(random.sample(bns, i))
    else:
        BNS = None
    BNS = COMM.bcast(BNS, root=MASTER_RANK)

    FC = FormCaller(mesh, space)
    f0 = FC('0-f', hybrid=ISH)
    f1 = FC('1-f', hybrid=ISH)

    f0.BC.boundaries = BNS[0]

    M1 = f1.matrices.mass
    E10 = f0.matrices.incidence
    BMAT = [[M1, E10], [E10.T, None]]
    SYSTEM = bmat(BMAT)
    SYSTEM.gathering_matrices = ([f1, f0], [f1, f0])

    SYSTEM.customize.identify_global_rows_according_to(1, f0.BC.interpret)

    SYSTEM = SYSTEM.assembled
    M = SYSTEM.do.gather_M_to_core()

    Gbd = f0.numbering.global_boundary_dofs

    if RANK == MASTER_RANK:
        bd = list()
        for bn in BNS[0]:
            bd.extend(Gbd[bn])

        bd = np.array(bd) + f1.num.global_dofs

        for i in range(f1.num.global_dofs + f0.num.global_dofs):
            if i in bd:
                assert M[i].nnz == 1, f"M[i].nnz={M[i].nnz}"
            else:
                assert M[i].nnz > 1

    # ------------ test for f2 -------------------------------------------------
    mesh, space = random_mesh_and_space_of_total_load_around(int(LOAD * 2 / 3),
                                                             exclude_periodic=True)
    # ------- test partial_dofs only having boundary dofs involved===
    bns = mesh.boundaries.names
    if RANK == MASTER_RANK:
        BNS = list()
        for _ in range(4):
            i = random.randint(1, len(bns))  # how many boundaries to be included.
            BNS.append(random.sample(bns, i))
    else:
        BNS = None
    BNS = COMM.bcast(BNS, root=MASTER_RANK)

    FC = FormCaller(mesh, space)
    f2 = FC('2-f', hybrid=ISH)
    f3 = FC('3-f', hybrid=ISH)

    f2.BC.boundaries = BNS[2]
    pd2 = f2.BC.interpret

    M3 = f3.matrices.mass
    E32 = f2.matrices.incidence
    BMAT = [[M3, E32], [E32.T, None]]
    SYSTEM = bmat(BMAT)
    SYSTEM.gathering_matrices = ([f3, f2], [f3, f2])
    SYSTEM.customize.identify_global_rows_according_to(1, pd2)
    SYSTEM = SYSTEM.assembled
    M = SYSTEM.do.gather_M_to_core()

    Gbd_f2 = f2.numbering.global_boundary_dofs

    if RANK == MASTER_RANK:
        bdf2 = list()
        for bn in BNS[2]:
            bdf2.extend(Gbd_f2[bn])

        bdf2 = np.array(bdf2) + f3.num.global_dofs

        for i in bdf2:
            assert M[i].nnz == 1, f"M[i].nnz={M[i].nnz}"
            assert M[i, i] == 1

    t2 = FC('2-t')
    t2.BC.boundaries = BNS[2]
    pdt = t2.BC.interpret

    M2 = f2.matrices.mass
    E32 = f2.matrices.incidence
    T2 = t2.matrices.trace
    BMAT = [[M2, E32.T, T2.T],
            [E32, None, None],
            [T2, None, None]]
    SYSTEM = bmat(BMAT)
    SYSTEM.gathering_matrices = ([f2, f3, t2], [f2, f3, t2])

    SYSTEM.customize.off_diagonally_identify_rows_according_to(
        2, 0, pdt, pd2)
    SYSTEM = SYSTEM.assembled
    M = SYSTEM.do.gather_M_to_core()

    GM_t2 = t2.numbering.gathering
    GM_F2 = f2.numbering.gathering

    gn_T2 = dict()
    gn_F2 = dict()
    for e in pdt.local.dofs:
        gn_T2[e] = GM_t2[e][pdt.local.dofs[e]]
        gn_F2[e] = GM_F2[e][pd2.local.dofs[e]]

    gn_T2 = COMM.gather(gn_T2, root=MASTER_RANK)
    gn_F2 = COMM.gather(gn_F2, root=MASTER_RANK)

    if RANK == MASTER_RANK:
        AT2 = dict()
        for _ in gn_T2:
            AT2.update(_)
        AF2 = dict()
        for _ in gn_F2:
            AF2.update(_)

        for e in AT2:
            gt = AT2[e]
            gf = AF2[e]

            gt = gt + f2.num.global_dofs + f3.num.global_dofs

            for i, j in zip(gt, gf):
                assert M[i].nnz == 1 and M[i, j] == 1

    # --------- single EWC -----------------------------------------------------------
    BMAT = [[M2, ], ]
    SYSTEM = bmat(BMAT)
    SYSTEM.gathering_matrices = ([f2, ], [f2, ])
    SYSTEM.customize.off_diagonally_identify_rows_according_to(
        0, 0, pd2, pd2)
    SYSTEM = SYSTEM.assembled
    M = SYSTEM.do.gather_M_to_core()

    if RANK == MASTER_RANK:
        # noinspection PyUnboundLocalVariable
        for e in AF2:
            gf = AF2[e]
            for i in gf:
                assert M[i].nnz == 1 and M[i, i] == 1
    # below, we check M2 is not changed by above adjusting.
    f2.BC.boundaries = BNS[3]
    pd2_new = f2.BC.interpret
    SYSTEM = bmat(BMAT)
    SYSTEM.gathering_matrices = ([f2, ], [f2, ])
    SYSTEM.customize.off_diagonally_identify_rows_according_to(
        0, 0, pd2_new, pd2_new)
    SYSTEM = SYSTEM.assembled
    M_new = SYSTEM.do.gather_M_to_core()
    gn_F2 = dict()
    for e in pd2_new.local.dofs:
        gn_F2[e] = GM_F2[e][pd2_new.local.dofs[e]]
    gn_F2 = COMM.gather(gn_F2, root=MASTER_RANK)
    if RANK == MASTER_RANK:
        AF2_new = dict()
        for _ in gn_F2:
            AF2_new.update(_)
        ALL = list()
        for e in AF2_new:
            ALL.extend(AF2_new[e])
        for i in range(np.shape(M_new)[0]):
            if i in ALL:
                assert M_new[i].nnz == 1 and M_new[i, i] == 1

    return 1


def test_TOOLS_NO14_partial_cochain_with_3dCSCG_form_BC():
    """"""
    if RANK == MASTER_RANK:
        print("-C- [test_TOOLS_NO14_partial_cochain_with_3dCSCG_form_BC] ....", flush=True)

    def Pressure(t, x, y, z): return 2.1 + t + np.cos(np.pi * x) * np.cos(2 * np.pi * y) * np.cos(3 * np.pi * z)

    def velocity_x(t, x, y, z):
        return 2.1 + t + np.cos(1.5*np.pi * x) * np.cos(2.5 * np.pi * y) * np.cos(3.5 * np.pi * z)

    def velocity_y(t, x, y, z):
        return 2.1 + t + np.cos(0.5*np.pi * x) * np.cos(2.1 * np.pi * y) * np.cos(3 * np.pi * z)

    def velocity_z(t, x, y, z):
        return 2.1 + t + np.cos(0.8*np.pi * x) * np.cos(1.3 * np.pi * y) * np.cos(0.7 * np.pi * z)

    if RANK == MASTER_RANK:
        ISH = (True, False, False, False)[random.randint(0, 3)]
        LOAD = random.randint(50, 85)
        time = random.random()
    else:
        ISH = None
        LOAD = None
        time = None
    ISH, LOAD, time = COMM.bcast([ISH, LOAD, time], root=MASTER_RANK)

    mesh, space = random_mesh_and_space_of_total_load_around(LOAD, exclude_periodic=True)

    bns = mesh.boundaries.names
    if RANK == MASTER_RANK:
        i = random.randint(1, len(bns))  # how many boundaries to be included.
        BNS = random.sample(bns, i)
    else:
        BNS = None
    BNS = COMM.bcast(BNS, root=MASTER_RANK)
    bcDs = dict()
    for bn in BNS:
        bcDs[bn] = Pressure
    bcDv = dict()
    for bn in BNS:
        bcDv[bn] = [velocity_x, velocity_y, velocity_z]

    FC = FormCaller(mesh, space)

    BS = FC('scalar', bcDs)
    SS = FC('scalar', Pressure)

    BV = FC('vector', bcDv)
    SV = FC('vector', [velocity_x, velocity_y, velocity_z])

    # ---  with 3d CSCG 0-form ---------------------------------------------------------------------
    f0 = FC('0-f', hybrid=ISH)
    f0.BC.CF = BS
    BS.current_time = time
    f0.BC.boundaries = BNS
    f0pc = f0.BC.interpret.local
    xi_et_sg = np.meshgrid(*space.nodes, indexing='ij')
    for i in f0pc.dofs:
        element = mesh.elements[i]
        local_dofs = f0pc.dofs[i]
        local_cochain = f0pc.cochains[i]
        x, y, z = element.coordinate_transformation.mapping(*xi_et_sg)
        x = x.ravel('F')[local_dofs]
        y = y.ravel('F')[local_dofs]
        z = z.ravel('F')[local_dofs]
        v_exact = Pressure(time, x, y, z)
        np.testing.assert_array_almost_equal(local_cochain, v_exact)
    #
    f0.BC.CF = SS
    SS.current_time = time
    f0.BC.boundaries = BNS
    f0pc = f0.BC.interpret.local
    xi_et_sg = np.meshgrid(*space.nodes, indexing='ij')
    for i in f0pc.dofs:
        element = mesh.elements[i]
        local_dofs = f0pc.dofs[i]
        local_cochain = f0pc.cochains[i]
        x, y, z = element.coordinate_transformation.mapping(*xi_et_sg)
        x = x.ravel('F')[local_dofs]
        y = y.ravel('F')[local_dofs]
        z = z.ravel('F')[local_dofs]
        v_exact = Pressure(time, x, y, z)
        np.testing.assert_array_almost_equal(local_cochain, v_exact)

    # ---  with 3d CSCG 2-form ---------------------------------------------------------------------
    f2 = FC('2-f', hybrid=ISH)
    f2.CF = SV
    SV.current_time = time
    f2.discretize()
    f2_cochain = f2.cochain.local

    f2.BC.CF = BV
    BV.current_time = time
    f2.BC.boundaries = BNS
    f2pc = f2.BC.interpret.local
    for i in f2pc.dofs:
        local_dofs = f2pc.dofs[i]
        local_cochain = f2pc.cochains[i]
        cochain_exact = f2_cochain[i][local_dofs]
        np.testing.assert_array_almost_equal(local_cochain, cochain_exact)

    f2.BC.CF = SV
    SV.current_time = time
    f2.BC.boundaries = BNS
    f2pc = f2.BC.interpret.local
    for i in f2pc.dofs:
        local_dofs = f2pc.dofs[i]
        local_cochain = f2pc.cochains[i]
        cochain_exact = f2_cochain[i][local_dofs]
        np.testing.assert_array_almost_equal(local_cochain, cochain_exact)

    # ---- with 3d CSCG 2-trace-form -----------------------------------------------
    t2 = FC('2-t')
    t2.CF = SV
    t2.discretize()
    t2_cochain = t2.cochain.local
    # 1: standard vector
    t2.BC.CF = SV
    t2.BC.boundaries = BNS
    t2pc = t2.BC.interpret.local
    for i in t2pc.dofs:
        local_dofs = t2pc.dofs[i]
        local_cochain = t2pc.cochains[i]
        cochain_exact = t2_cochain[i][local_dofs]
        np.testing.assert_array_almost_equal(local_cochain, cochain_exact)
    # 2: boundary-wise vector
    t2.BC.CF = BV
    t2.BC.boundaries = BNS
    t2pc = t2.BC.interpret.local
    for i in t2pc.dofs:
        local_dofs = t2pc.dofs[i]
        local_cochain = t2pc.cochains[i]
        cochain_exact = t2_cochain[i][local_dofs]
        np.testing.assert_array_almost_equal(local_cochain, cochain_exact)

    t2.CF = SS
    t2.discretize()
    t2_cochain = t2.cochain.local
    # 1: standard scalar
    t2.BC.CF = SS
    t2.BC.boundaries = BNS
    t2pc = t2.BC.interpret.local
    for i in t2pc.dofs:
        local_dofs = t2pc.dofs[i]
        local_cochain = t2pc.cochains[i]
        cochain_exact = t2_cochain[i][local_dofs]
        np.testing.assert_array_almost_equal(local_cochain, cochain_exact)
    # 2: boundary-wise scalar
    t2.BC.CF = BS
    t2.BC.boundaries = BNS
    t2pc = t2.BC.interpret.local
    for i in t2pc.dofs:
        local_dofs = t2pc.dofs[i]
        local_cochain = t2pc.cochains[i]
        cochain_exact = t2_cochain[i][local_dofs]
        np.testing.assert_array_almost_equal(local_cochain, cochain_exact)

    # test set_entries_according_to_CSCG_partial_cochains for EWC vectors.
    cf0 = EWC_ColumnVector(mesh, f0.num.basis)
    cf2 = EWC_ColumnVector(mesh, f2.num.basis)
    ct2 = EWC_ColumnVector(mesh, t2.num.basis)

    t2pc = t2.BC.interpret
    b = concatenate([cf0, cf2, ct2])
    b.gathering_matrix = [f0, f2, t2]
    b.customize.set_entries_according_to(2, t2pc)
    B = b.assembled
    B = B.do.gather_V_to_core()

    GM = t2.numbering.gathering
    gn_t2 = dict()
    for e in t2pc.local.dofs:
        gn_t2[e] = GM[e][t2pc.local.dofs[e]]
    gn_t2 = COMM.gather(gn_t2, root=MASTER_RANK)
    if RANK == MASTER_RANK:
        At2 = dict()
        for _ in gn_t2:
            At2.update(_)

        ALL = list()
        for e in At2:
            ALL.extend(At2[e])
        ALL = np.array(ALL) + f0.num.global_dofs + f2.num.global_dofs

        for i, bi in enumerate(B):
            if i in ALL:
                assert bi != 0
            else:
                assert bi == 0

    b.customize.set_constant_entries_according_to(2, t2pc, -1)
    B = b.assembled
    B = B.do.gather_V_to_core()
    if RANK == MASTER_RANK:
        for i, bi in enumerate(B):
            # noinspection PyUnboundLocalVariable
            if i in ALL:
                assert bi == -1
            else:
                assert bi == 0

    b = concatenate([ct2, cf2, cf0])
    b.gathering_matrix = [t2, f2, f0]
    b.customize.set_constant_entries_according_to(0, t2pc, -2.1415)
    B = b.assembled
    B = B.do.gather_V_to_core()
    if RANK == MASTER_RANK:
        ALL = ALL - f0.num.global_dofs - f2.num.global_dofs
        for i, bi in enumerate(B):
            if i in ALL:
                assert bi == -2.1415
            else:
                assert bi == 0

    return 1


def test_TOOLS_NO15_linear_system_apply_BC():
    """"""
    if RANK == MASTER_RANK:
        LOAD = random.randint(50, 300)
        time = random.random()
        print(f"-S- [test_TOOLS_NO15_linear_system_apply_BC] @ load = {LOAD}, time=%.2f..." % time, flush=True)
    else:
        LOAD = None
        time = None
    LOAD, time = COMM.bcast([LOAD, time], root=MASTER_RANK)

    mesh, space = random_mesh_and_space_of_total_load_around(LOAD, exclude_periodic=True, mesh_boundary_num='>=2')
    FC = FormCaller(mesh, space)

    def Pressure(t, x, y, z):
        return 2.5 + t + np.cos(np.pi * x) * np.cos(2 * np.pi * y) * np.cos(3 * np.pi * z)

    def velocity_x(t, x, y, z):
        return 2.5 + t + np.cos(1.5*np.pi * x) * np.cos(2.5 * np.pi * y) * np.cos(3.5 * np.pi * z)

    def velocity_y(t, x, y, z):
        return 2.5 + t + np.cos(0.5*np.pi * x) * np.cos(2.1 * np.pi * y) * np.cos(3.2 * np.pi * z)

    def velocity_z(t, x, y, z):
        return 2.5 + t + np.cos(0.8*np.pi * x) * np.cos(1.3 * np.pi * y) * np.cos(0.7 * np.pi * z)
    # SS = FC('scalar', Pressure)
    # SV = FC('vector', [velocity_x, velocity_y, velocity_z])

    bns = mesh.boundaries.names
    if RANK == MASTER_RANK:
        i = random.randint(1, len(bns)-1)  # how many boundaries to be included.
        BNS = random.sample(bns, i)
    else:
        BNS = None
    BNS = COMM.bcast(BNS, root=MASTER_RANK)
    bcDs = dict()
    for bn in BNS:
        bcDs[bn] = Pressure
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
    for bn in BNS_com:
        bcDv_com[bn] = [velocity_x, velocity_y, velocity_z]
    # BS_com = FC('scalar', bcDs_com)
    BV_com = FC('vector', bcDv_com)

    # Poisson hybrid system ---------------------------------------------------------------------
    f2 = FC('2-f', hybrid=True)
    f3 = FC('3-f', hybrid=True)
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

    BS.current_time = time
    t2.BC.CF = BS
    t2.BC.boundaries = BNS
    t2BC1 = t2.BC.interpret

    Axb.customize.apply_strong_BC(2, 2, t2BC1)

    aA, ab = Axb.assembled
    aA = aA.do.gather_M_to_core()
    ab = ab.do.gather_V_to_core()

    GM_t2 = t2.numbering.gathering

    t2GBD = t2.numbering.global_boundary_dofs
    if RANK == MASTER_RANK:

        dofs_changed = list()  # the rows that been changed in the global matrix.
        for bn in BNS:
            dofs_changed.extend(t2GBD[bn])
        dofs_changed = np.array(dofs_changed) + f2.num.global_dofs + f3.num.global_dofs
        for i in dofs_changed:
            assert aA[i].nnz == 1 and aA[i, i] == 1 and ab[i] != 0

        NOT_changed = list()  # the rows that not been changed in the global matrix.
        for i in range(aA.shape[0]):
            if i in dofs_changed:
                pass
            else:
                NOT_changed.append(i)

        ab_SUB = ab[NOT_changed]
        assert np.all(ab_SUB == 0), f"not change places must be all zero!"
        aA_SUB = aA[NOT_changed, :]
        aA_SUB = aA_SUB[f2.num.global_dofs:, f2.num.global_dofs:]
        assert aA_SUB.nnz == 0, f"not change places must be all zero!"

    gn_T2 = dict()
    for e in t2BC1.local.dofs:
        gn_T2[e] = GM_t2[e][t2BC1.local.dofs[e]]

    gn_T2 = COMM.gather(gn_T2, root=MASTER_RANK)

    if RANK == MASTER_RANK:
        AT2 = dict()
        for _ in gn_T2:
            AT2.update(_)

        changed = list()

        for e in AT2:
            gt = AT2[e]
            gt = gt + f2.num.global_dofs + f3.num.global_dofs

            for i in gt:
                assert aA[i].nnz == 1 and aA[i, i] == 1 and ab[i] != 0

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
        assert np.all(ab == 0), f"not change places must be all zero!"
        aA = aA[not_changed, :]
        aA = aA[f2.num.global_dofs:, f2.num.global_dofs:]
        assert aA.nnz == 0, f"not change places must be all zero!"

    t2.BC.CF = None
    t2.BC.boundaries = BNS_com
    t2BC2 = t2.BC.interpret
    BV_com.current_time = time
    f2.BC.CF = BV_com
    f2.BC.boundaries = BNS_com
    f2BC2 = f2.BC.interpret

    Axb.customize.apply_strong_BC(2, 0, t2BC2, f2BC2)

    aA, ab = Axb.assembled
    aA = aA.do.gather_M_to_core()
    ab = ab.do.gather_V_to_core()

    if RANK == MASTER_RANK:

        dofs_changed = list()  # the rows that been changed in the global matrix.
        for bn in t2GBD:  # all the boundary dofs.
            dofs_changed.extend(t2GBD[bn])
        dofs_changed = np.array(dofs_changed) + f2.num.global_dofs + f3.num.global_dofs

        dofs_changed_2 = list()  # the rows that been changed in the global matrix in the second apply_strong_BC
        for bn in BNS_com:  # all the boundary dofs.
            dofs_changed_2.extend(t2GBD[bn])
        dofs_changed_2 = np.array(dofs_changed_2) + f2.num.global_dofs + f3.num.global_dofs

        for i in dofs_changed:
            assert aA[i].nnz == 1 and ab[i] != 0

            if i in dofs_changed_2:
                assert aA[i, i] == 0

        NOT_changed = list()  # the rows that not been changed in the global matrix.
        for i in range(aA.shape[0]):
            if i in dofs_changed:
                pass
            else:
                NOT_changed.append(i)
        ab_SUB = ab[NOT_changed]
        assert np.all(ab_SUB == 0), f"not change places must be all zero!"
        aA_SUB = aA[NOT_changed, :]
        aA_SUB = aA_SUB[f2.num.global_dofs:, f2.num.global_dofs:]
        assert aA_SUB.nnz == 0, f"not change places must be all zero!"

        # noinspection PyUnboundLocalVariable
        for i in changed:  # the changed rows after the first time `apply_strong_BC`.
            assert aA[i].nnz == 1 and aA[i, i] == 1 and ab[i] != 0

    # check changes in block[2][0] .............................
    GM_F2 = f2.numbering.gathering
    gn_T2 = dict()
    gn_F2 = dict()
    for e in t2BC2.local.dofs:
        gn_T2[e] = GM_t2[e][t2BC2.local.dofs[e]]
        gn_F2[e] = GM_F2[e][f2BC2.local.dofs[e]]
    gn_T2 = COMM.gather(gn_T2, root=MASTER_RANK)
    gn_F2 = COMM.gather(gn_F2, root=MASTER_RANK)
    if RANK == MASTER_RANK:
        AT2 = dict()
        for _ in gn_T2:
            AT2.update(_)
        AF2 = dict()
        for _ in gn_F2:
            AF2.update(_)
        for e in AT2:
            gt = AT2[e]
            gf = AF2[e]
            gt = gt + f2.num.global_dofs + f3.num.global_dofs
            for i, j in zip(gt, gf):
                assert aA[i].nnz == 1 and aA[i, j] == 1 and ab[i] != 0
                assert i not in NOT_changed, f"dof #{i} must be changed"

    return 1


if __name__ == '__main__':
    # mpiexec -n 5 python
    pass
