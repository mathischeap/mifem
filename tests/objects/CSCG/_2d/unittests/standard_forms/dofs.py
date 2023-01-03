# -*- coding: utf-8 -*-
"""test dofs related functions."""
import sys
if './' not in sys.path:
    sys.path.append('./')
from root.config.main import *
from tests.objects.CSCG._2d.randObj.form_caller import random_FormCaller_of_total_load_around
import random


def test_standard_forms_DOFS():
    """"""
    if RANK == MASTER_RANK:
        load = random.randint(10, 100)
        IH = [True, False][random.randint(0, 1)]
        print(f"~~~ [test_standard_forms_DOFS] @ load= {load}... ", flush=True)
    else:
        load = None
        IH = None
    load, IH = COMM.bcast([load, IH], root=MASTER_RANK)

    FC = random_FormCaller_of_total_load_around(load, exclude_periodic=True)

    region_name = FC._mesh_.domain.regions.names[0]
    def P(t, x, y): return - np.sin(0.715*np.pi*x) * np.cos(0.8451*np.pi*y) + t/3
    def Q(t, x, y): return np.cos(1.1215*np.pi*x) * np.sin(0.121*np.pi*y) - t/3
    scalar = FC('scalar', P)
    vector = FC('vector', (P, Q))

    F0 = [FC('0-f-o', hybrid=IH), FC('0-f-i', hybrid=IH)]
    F1 = [FC('1-f-o', hybrid=IH), FC('1-f-i', hybrid=IH)]
    F2 = [FC('2-f-o', hybrid=IH), FC('2-f-i', hybrid=IH)]

    for f in F0:
        dof = f.dofs.do.find.dof_at_corner_of_region(region_name, 'RD')

    return 1




if __name__ == '__main__':
    # mpiexec -n 4 python objects\CSCG\_2d\__tests__\unittests\standard_forms\dofs.py
    test_standard_forms_DOFS()