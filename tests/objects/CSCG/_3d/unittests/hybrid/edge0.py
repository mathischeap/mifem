# -*- coding: utf-8 -*-



import sys
if './' not in sys.path: sys.path.append('./')
from root.config.main import RANK, MASTER_RANK, COMM
import random
from tools.linearAlgebra.elementwiseCache.operators.bmat.main import bmat
from tools.linearAlgebra.elementwiseCache.objects.sparseMatrix.main import EWC_SparseMatrix
from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller





def test_hybridization_of_standard_0_form():
    """"""
    if RANK == MASTER_RANK:
        mesh_no = random.randint(0,1)
        print(f"H-0 [test_hybridization_of_standard_0_form] @mesh_no={mesh_no}... ", flush=True)
    else:
        mesh_no = None

    mesh_no = COMM.bcast(mesh_no, root=MASTER_RANK)

    if mesh_no == 0:
        mesh = MeshGenerator('crazy', c=0.1)([2,3,4])
    elif mesh_no == 1:
        mesh = MeshGenerator('cuboid', region_layout=[3, 2, 1])([1,2,3])
    else:
        raise NotImplementedError(f"mesh_no={mesh_no} not implemented.")

    bns = mesh.boundaries.names
    if RANK == MASTER_RANK:
        Dirichlet_boundaries = random.sample(bns, random.randint(1,5))
    else:
        Dirichlet_boundaries = None
    Dirichlet_boundaries = COMM.bcast(Dirichlet_boundaries, root=MASTER_RANK)
    Neumann_boundaries = list()
    for bn in bns:
        if bn not in Dirichlet_boundaries:
            Neumann_boundaries.append(bn)
    space = SpaceInvoker('polynomials')([2, 1, 3])
    def pressure(t, x, y, z): return 100 + 0 * x * y * z * t
    FC = FormCaller(mesh, space)
    pressure = FC('scalar', pressure)

    f = FC('0-f', is_hybrid=True)
    t = FC('0-adt')
    e = FC('0-e')

    f.BC.CF = pressure
    f.BC.CF.current_time = 0
    f.BC.boundaries = Neumann_boundaries

    t.prime.BC.CF = pressure
    t.prime.BC.CF.current_time = 0
    t.BC.boundaries = Dirichlet_boundaries

    T, D, C, b, eGM = f.special.hybrid_pairing(t, e)

    for i in range(mesh.node.elements.GLOBAL_num):
        f.dofs.visualize.matplot.connection_around_node_element(i, T, D, C, t, e, checking_mode=True)

    Id = EWC_SparseMatrix(mesh, ('identity', f.num.basis))
    T_T = t.matrices.trace.T
    A = ([Id,   T_T, None],
         [ T,     D,    C],
         [None, C.T, None])
    A = bmat(A)
    A.gathering_matrices = ([f, t, eGM], [f, t, eGM])
    A = A.assembled
    assert A.condition.condition_number < 100, f"We should get a non-singular system."




    return 1





if __name__ == '__main__':
    # mpiexec -n 4 python objects\CSCG\_3d\__tests__\unittests\hybrid\edge0.py
    test_hybridization_of_standard_0_form()