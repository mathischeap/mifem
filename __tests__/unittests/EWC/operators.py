

import sys
if './' not in sys.path: sys.path.append('./')
from objects.CSCG._2d.__tests__.Random.form_caller import random_FormCaller_of_total_load_around as _2d_RANDOM_FC_
import random

from numpy import sin, cos, exp, pi

from root.config.main import rAnk, mAster_rank, np, cOmm



def test_LinearAlgebra_EWC_No1_Operators():
    """"""
    #--------- 2d CSCG tests ---------------------------------------------------------------
    if rAnk == mAster_rank:
        print("&&& [test_LinearAlgebra_EWC_No1_Operators] ...... ", flush=True)
        load = random.randint(10, 999)
        IH = [True, False][random.randint(0, 1)]
    else:
        load = None
        IH = None
    load, IH = cOmm.bcast([load, IH], root=mAster_rank)

    # FC = random_FormCaller_of_total_load_around(load, mesh_pool='rectangle_periodic')
    FC = _2d_RANDOM_FC_(load)
    nu = random.random()

    a, b, c, d = random.random(), random.random(), random.random(), random.random()
    rT = random.random()
    def P(t, x, y): return - sin(a*pi*x) * cos(b*pi*y) * exp(- 0.2 * pi * t )
    # def Q(t, x, y): return cos(c*pi*x) * sin(d*pi*y) * exp(- 0.3 * pi * t )

    scalar = FC('scalar', P)

    w = FC('0-f-o', is_hybrid=False, name='vorticity')
    u = FC('1-f-o', is_hybrid=False, name='velocity')

    M1 = u.matrices.mass
    E10 = w.matrices.incidence

    w.TW.func.do.set_func_body_as(scalar)
    w.TW.current_time = rT
    w.TW.do.push_all_to_instant()
    w.discretize()

    A = nu * M1 @ E10 @ w
    B = nu * (M1 @ E10) @ w

    for i in A:
        ai = nu * M1[i].toarray() @ (E10[i].toarray() @ w.cochain.EWC[i].toarray())
        Ai = A[i].toarray()
        Bi = B[i].toarray()
        np.testing.assert_array_almost_equal(Ai - ai, 0, decimal=5)
        np.testing.assert_array_almost_equal(Bi - ai, 0, decimal=5)



    return 1

if __name__ == '__main__':
    # mpiexec -n 4 python __tests__\unittests\EWC\operators.py

    test_LinearAlgebra_EWC_No1_Operators()
