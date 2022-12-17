# -*- coding: utf-8 -*-

import sys
if './' not in sys.path:
    sys.path.append('./')
from tests.objects.CSCG._2d.randObj.form_caller import random_FormCaller_of_total_load_around as _2d_RANDOM_FC_
import random

from numpy import sin, cos, exp, pi

from root.config.main import RANK, MASTER_RANK, np, COMM


def test_LinearAlgebra_EWC_No1_Operators():
    """"""
    # -------- 2d CSCG tests ---------------------------------------------------------------
    if RANK == MASTER_RANK:
        print("&&& [test_LinearAlgebra_EWC_No1_Operators] ...... ", flush=True)
        load = random.randint(180, 199)
        IH = [True, False][random.randint(0, 1)]
    else:
        load = None
        IH = None
    load, IH = COMM.bcast([load, IH], root=MASTER_RANK)

    # FC = random_FormCaller_of_total_load_around(load, mesh_pool='rectangle_periodic')
    FC = _2d_RANDOM_FC_(load)
    nu = random.random()

    a, b, c, d = random.random(), random.random(), random.random(), random.random()
    rT = random.random()
    def P(t, x, y): return - sin(a*pi*x) * cos(b*pi*y) * exp(- 0.2 * pi * t)

    scalar = FC('scalar', P)

    w = FC('0-f-o', hybrid=False, name='vorticity')
    u = FC('1-f-o', hybrid=False, name='velocity')

    M1 = u.matrices.mass
    E10 = w.matrices.incidence

    w.CF = scalar
    w.CF.current_time = rT
    w.discretize()

    A = nu * M1 @ E10 @ w
    B = nu * (M1 @ E10) @ w

    for i in A:
        ai = nu * M1[i].toarray() @ (E10[i].toarray() @ w.cochain.EWC[i].toarray())
        Ai = A[i].toarray()
        Bi = B[i].toarray()
        np.testing.assert_array_almost_equal(Ai, ai, decimal=3)
        np.testing.assert_array_almost_equal(Bi, ai, decimal=3)

    return 1


if __name__ == '__main__':
    # mpiexec -n 4 python tests/tools/EWC/operators.py

    test_LinearAlgebra_EWC_No1_Operators()
