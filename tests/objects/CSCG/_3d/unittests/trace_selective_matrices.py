# -*- coding: utf-8 -*-
import sys
if './' not in sys.path:
    sys.path.append('./')

from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller
import random
from root.config.main import RANK, MASTER_RANK, np, COMM
from scipy.sparse import csc_matrix


def test_trace_and_selective_matrices():
    """"""
    if RANK == MASTER_RANK:
        print(f"T~S [test_trace_and_selective_matrices] ...... ", flush=True)
        i, j, k = random.randint(1, 4), random.randint(1, 4), random.randint(1, 4)
    else:
        i, j, k = None, None, None
    i, j, k = COMM.bcast([i, j, k], root=MASTER_RANK)

    mesh = MeshGenerator('crazy')([k, j, i])
    space = SpaceInvoker('polynomials')([('Lobatto', i), ('Lobatto', j), ('Lobatto', k)])
    FC = FormCaller(mesh, space)

    f0 = FC('0-f')
    f1 = FC('1-f')
    f2 = FC('2-f')

    t0 = FC('0-t')
    t1 = FC('1-t')
    t2 = FC('2-t')

    S0 = t0.matrices.selective
    S1 = t1.matrices.selective
    S2 = t2.matrices.selective

    T0 = t0.matrices.trace
    T1 = t1.matrices.trace
    T2 = t2.matrices.trace

    s0f = list()
    for side in 'NSWEBF':
        s0f.append(f0.numbering.do.find.local_dofs_on_element_side(side))
    s0f = np.concatenate(s0f)
    s1f = list()
    for side in 'NSWEBF':
        s1f.append(f1.numbering.do.find.local_dofs_on_element_side(side))
    s1f = np.concatenate(s1f)
    s2f = list()
    for side in 'NSWEBF':
        s2f.append(f2.numbering.do.find.local_dofs_on_element_side(side))
    s2f = np.concatenate(s2f)

    SS0 = csc_matrix((np.ones(t0.num.basis), ([_ for _ in range(t0.num.basis)], s0f)),
                     shape=(t0.num.basis, f0.num.basis))
    SS1 = csc_matrix((np.ones(t1.num.basis), ([_ for _ in range(t1.num.basis)], s1f)),
                     shape=(t1.num.basis, f1.num.basis))
    SS2 = csc_matrix((np.ones(t2.num.basis), ([_ for _ in range(t2.num.basis)], s2f)),
                     shape=(t2.num.basis, f2.num.basis))
    for e in S0:
        assert (S0[e] - SS0).nnz == 0
        assert (S1[e] - SS1).nnz == 0
        assert (S2[e] - SS2).nnz == 0
        break  # only test once as same in all mesh elements.

    for e in T0:
        t = T0[e].tocsc()
        assert np.all(t.indices == SS0.indices)
        assert np.all(t.indptr == SS0.indptr)
        t = T1[e].tocsc()
        assert np.all(t.indices == SS1.indices)
        assert np.all(t.indptr == SS1.indptr)
        t = T2[e].tocsc()
        assert np.all(t.indices == SS2.indices)
        assert np.all(t.indptr == SS2.indptr)
        break  # only test once as same in all mesh elements.

    return 1


if __name__ == '__main__':
    # mpiexec -n 4 python objects\CSCG\_3d\__tests__\unittests\trace_selective_matrices.py
    test_trace_and_selective_matrices()
