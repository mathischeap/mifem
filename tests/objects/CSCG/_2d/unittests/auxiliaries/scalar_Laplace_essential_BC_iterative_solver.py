# -*- coding: utf-8 -*-

from tools.miLinearAlgebra.linearSystem.main import LinearSystem
from objects.CSCG._2d.master import MeshGenerator, SpaceInvoker, FormCaller, ExactSolutionSelector
from tools.elementwiseCache.dataStructures.objects.columnVector.main import EWC_ColumnVector
from tools.elementwiseCache.dataStructures.operators.bmat.main import bmat
from tools.elementwiseCache.dataStructures.operators.concatenate.main import concatenate


def scalar_Laplace_solver_iterative_solver(c, Kx, Ky, Nx, Ny):
    """

    Parameters
    ----------
    c
    Kx
    Ky
    Nx
    Ny

    Returns
    -------

    """
    mesh = MeshGenerator('crazy', c=c, bounds=[(0,1), (0,1)])([Kx, Ky], EDM=None)
    space = SpaceInvoker('polynomials')([('Lobatto', Nx), ('Lobatto', Ny)])
    FC = FormCaller(mesh, space)
    ES = ExactSolutionSelector(mesh)('sL:sincos1')

    u = FC('1-f-o', hybrid=False)
    p = FC('2-f-o', hybrid=False)
    f = FC('2-f-o', hybrid=False)

    B0 = EWC_ColumnVector(mesh, u.num.basis)
    B0.gathering_matrix = u

    f.CF = ES.source_term
    f.CF.current_time = 0
    f.discretize()
    B1 = f.cochain.EWC
    B1.gathering_matrix = p

    M1 = u.matrices.mass
    M2 = p.matrices.mass
    E21 = u.matrices.incidence
    E12 = E21.T
    E12M2 = E12 @ M2

    A = bmat(([M1 , E12M2],
              [E21, None ]))
    A.gathering_matrices = ((u, p), (u, p))
    b = concatenate([B0, B1])
    LS = LinearSystem(A, b)

    results = LS.solve('GMRES')(0, restart=100, maxiter="20")[0]
    results.do.distributed_to(u, p)

    u.CF = ES.velocity
    u.CF.current_time = 0
    u_error_L2 = u.error.L()
    p.CF = ES.potential
    p.CF.current_time = 0
    p_error_L2 = p.error.L()

    return u_error_L2, p_error_L2