# -*- coding: utf-8 -*-
import sys
if './' not in sys.path: sys.path.append('./')

from objects.CSCG._3d.__init__ import mesh as mesh3
from objects.CSCG._3d.__init__ import space as space3
from objects.CSCG._3d.__init__ import form as form3

from root.config.main import RANK, MASTER_RANK, COMM
import random

def test_hybridization_trace2():
    """"""
    mesh = mesh3('cuboid', region_layout=[2,2,2])([3,2,1])
    space = space3('polynomials')([('Lobatto', 1), ('Lobatto', 2), ('Lobatto', 3)])
    FC = form3(mesh, space)

    all_boundaries = mesh.boundaries.names
    if RANK == MASTER_RANK:
        rn = random.randint(1,5)
        boundaries = random.sample(all_boundaries, rn)
    else:
        boundaries = None
    boundaries = COMM.bcast(boundaries, root=MASTER_RANK)

    u_boundaries = boundaries
    if RANK == MASTER_RANK:
        print(f"HT2 [test_hybridization_trace2]...", flush=True)

    p_boundaries = list()
    for b in all_boundaries:
        if b not in u_boundaries:
            p_boundaries.append(b)

    def ux_func(t, x, y, z): return 100 + t + x + y + z
    def uy_func(t, x, y, z): return 100 + t + x + y + z
    def uz_func(t, x, y, z): return 100 + t + x + y + z
    velocity = FC('vector', [ux_func, uy_func, uz_func])

    def p_func(t, x, y, z): return 10 + t + x + y + z
    pressure = FC('scalar', p_func)

    u = FC('2-f', is_hybrid=True, name='velocity')
    t = FC('2-adt', name='pressure_trace')

    u.BC.CF = velocity
    u.BC.boundaries = u_boundaries
    t.prime.BC.CF = pressure
    t.BC.boundaries = p_boundaries

    T, D, b = u.special.hybrid_pairing(t, time=0)

    u.dofs.do.hybrid_pairing_check(t, T, D, b)

    return 1


if __name__ == '__main__':
    # mpiexec -n 4 python objects/CSCG/_3d/__tests__/unittests/hybrid/trace2.py

    test_hybridization_trace2()