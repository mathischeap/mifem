# -*- coding: utf-8 -*-

import sys
if './' not in sys.path: sys.path.append('./')
from root.config.main import rAnk, mAster_rank, cOmm
import random
from objects.CSCG._3d.__tests__.Random.form_caller import random_FormCaller_of_total_load_around
from screws.miscellaneous.mios import mkdir, remove, rmdir
from screws.miscellaneous.random_string.digits import randomStringDigits
from tools.linear_algebra.elementwise_cache.operators.bmat.main import bmat
from tools.linear_algebra.elementwise_cache.objects.sparse_matrix.main import EWC_SparseMatrix
from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller, ExactSolutionSelector




def test_hybridization_of_standard_1_form():
    """"""
    if rAnk == mAster_rank:
        load = random.randint(30,75)
        print(f"H-1 [test_hybridization_of_standard_1_form] @load={load}... ", flush=True)
        image_folder = '__hybrid_1_test_' + randomStringDigits(5) +'__'
    else:
        load = None
        image_folder = None
    load = cOmm.bcast(load, root=mAster_rank)
    image_folder = cOmm.bcast(image_folder, root=mAster_rank)
    FC = random_FormCaller_of_total_load_around(load)
    mesh = FC.mesh

    f1 = FC('1-f', is_hybrid=True)
    t1 = FC('1-adt')
    e1 = FC('1-e')

    #-------------- test 1 -------------------------------------------------------------------
    T = t1.matrices.trace
    C = e1.matrices.complement

    Id = EWC_SparseMatrix(mesh, ('identity', f1.num.basis))

    mkdir(image_folder)
    t1N = t1.prime.numbering.gathering.GLOBAL_num_dofs
    for i in range(t1N):
        f1.dofs.visualize.matplot.connection_through_trace_dof(
            i, T, C, t1, e1, checking_mode=True)

    if rAnk == mAster_rank:
        I = random.sample(range(t1N), 5)
    else:
        I = None
    I = cOmm.bcast(I, root=mAster_rank)
    for i in I:
        f1.dofs.visualize.matplot.connection_through_trace_dof(
            i, T, C, t1, e1, checking_mode=False,
            saveto=image_folder + f'/{i}.png')
    for i in I:
        remove(image_folder + f'/{i}.png')

    #-------------- test 2 -------------------------------------------------------------------
    T, C = f1.special.___PRIVATE_overcoming_hybrid_singularity___(T, C)[:2]
    e1N = e1.numbering.gathering.GLOBAL_num_dofs
    for i in range(e1N):
        f1.dofs.visualize.matplot.connection_through_around_edge_dof(
            i, T, C, t1, e1, checking_mode=True)
    if rAnk == mAster_rank:
        I = random.sample(range(e1N), 5)
    else:
        I = None
    I = cOmm.bcast(I, root=mAster_rank)
    for i in I:
        f1.dofs.visualize.matplot.connection_through_around_edge_dof(
            i, T, C, t1, e1, checking_mode=False,
            saveto=image_folder + f'/{i}.png')
    for i in I:
        remove(image_folder + f'/{i}.png')

    rmdir(image_folder)

    A = ([Id  ,  T.T, None],
         [T   , None, C   ],
         [None,  C.T, None])
    A = bmat(A)
    A.gathering_matrices = ([f1, t1, e1], [f1, t1, e1])
    A = A.assembled
    assert A.condition.condition_number < 1000, f"no longer singular!"

    #-------------- test 3 -------------------------------------------------------------------
    mesh = MeshGenerator('crazy', c=0.1)([2,3,4])
    ES = ExactSolutionSelector(mesh)('Poisson:sincos1')
    if rAnk == mAster_rank:
        Dirichlet_boundaries = random.sample(
            ['Back', 'Front', 'West', "East", 'South', 'North'], random.randint(0,6)) #
    else:
        Dirichlet_boundaries = None
    Dirichlet_boundaries = cOmm.bcast(Dirichlet_boundaries, root=mAster_rank)
    Neumann_boundaries = list()
    for bn in mesh.boundaries.names:
        if bn not in Dirichlet_boundaries:
            Neumann_boundaries.append(bn)

    space = SpaceInvoker('polynomials')([4, 3, 2])
    FC = FormCaller(mesh, space)

    f1 = FC('1-f', is_hybrid=True)
    t1 = FC('1-adt')
    e1 = FC('1-e')

    f1.TW.BC.body = ES.status.velocity
    f1.TW.do.push_BC_to_instant(0)
    f1.BC.valid_boundaries = Neumann_boundaries

    t1.prime.TW.BC.body = ES.status.velocity.components.T_perp
    t1.prime.TW.do.push_BC_to_instant(0)
    t1.BC.valid_boundaries = Dirichlet_boundaries

    T, D, C, b = f1.special.hybrid_pairing(t1, e1)[:4]
    e1N = e1.numbering.gathering.GLOBAL_num_dofs
    for i in range(e1N):
        f1.dofs.visualize.matplot.connection_through_around_edge_dof(
            i, T, C, t1, e1, checking_mode=True)

    return 1













if __name__ == '__main__':
    # mpiexec -n 4 python objects\CSCG\_3d\__tests__\unittests\hybrid\edge1.py
    test_hybridization_of_standard_1_form()