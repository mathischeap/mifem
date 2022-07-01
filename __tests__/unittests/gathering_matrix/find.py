# -*- coding: utf-8 -*-
"""Here we test the find function of Gathering matrix."""

import sys
if './' not in sys.path: sys.path.append('./')


from numpy import sin, cos, exp, pi
from root.config.main import *
import random
from objects.CSCG._3d.__tests__.Random.form_caller import random_FormCaller_of_total_load_around as _3d_Caller_
from objects.CSCG._2d.__tests__.Random.form_caller import random_FormCaller_of_total_load_around as _2d_Caller_



def test_GatheringMatrix_find():
    """"""
    if rAnk == mAster_rank:
        print("&F& [test_GatheringMatrix_find] ...... ", flush=True)

    if rAnk == mAster_rank:
        load = random.randint(30, 200)
        IH = [True, False][random.randint(0,1)]
    else:
        load = None
        IH = None
    load, IH = cOmm.bcast([load, IH], root=mAster_rank)

    if rAnk == mAster_rank:
        a, b, c, d = random.random(), random.random(), random.random(), random.random()
        rT = random.random()
    else:
        a, b, c, d, rT = [None for _ in range(5)]

    a, b, c, d, rT = cOmm.bcast([a, b, c, d, rT], root=mAster_rank)

    def func2(t, x, y):
        return - sin(a * pi * x) * cos(b * pi * y) * exp(d * pi * t / 10) + (a+b+c+d)/4
    def func3(t, x, y, z):
        return - sin(a * pi * x) * cos(b * pi * y) * sin( c * pi * z) * exp(d * pi * t / 10) + (a+b+c+d)/4

    RC = [_2d_Caller_, _3d_Caller_]
    for rc in RC:
        FC = rc(load, exclude_periodic=True)
        if FC._mesh_.ndim == 2:
            f = FC('0-f-o', is_hybrid=False, name='vorticity')
            scalar = FC('scalar', func2)
        else:
            f = FC('0-f', is_hybrid=False, name='vorticity')
            scalar = FC('scalar', func3)

        f.TW.func.do.set_func_body_as(scalar)
        f.TW.current_time = rT
        f.TW.do.push_all_to_instant()
        f.discretize()

        GLOBAL_num_dofs = f.num.GLOBAL_dofs
        if rAnk == mAster_rank:
            num_samples = random.randint(1, int(load/5))
            dofs = random.sample(range(GLOBAL_num_dofs), num_samples)
        else:
            dofs = None
        dofs = cOmm.bcast(dofs, root=mAster_rank)

        globe = f.cochain.globe
        for dof in dofs:
            dof_cc = f.cochain.dofwise[dof]
            if globe.V[dof,0] !=0:
                np.testing.assert_almost_equal(dof_cc, globe.V[dof,0])

    return 1




if __name__ == '__main__':
    # mpiexec -n 4 python tests\unittests\gathering_matrix\find.py

    test_GatheringMatrix_find()
