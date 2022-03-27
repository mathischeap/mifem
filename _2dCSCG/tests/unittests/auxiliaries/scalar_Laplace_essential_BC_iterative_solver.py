
# mpiexec -n 4 python _2dCSCG\tests\unittests\auxiliaries\scalar_Laplace_iterative_solver.py

import sys
if './' not in sys.path: sys.path.append('./')
from tools.linear_algebra.linear_system.main import LinearSystem
from _2dCSCG.master import MeshGenerator, SpaceInvoker, FormCaller, ExactSolutionSelector
from tools.linear_algebra.elementwise_cache.objects.column_vector.main import EWC_ColumnVector
from tools.linear_algebra.elementwise_cache.operators.bmat.main import bmat
from tools.linear_algebra.elementwise_cache.operators.concatenate.main import concatenate





def scalar_Laplace_solver_iterative_solver(c, Kx, Ky, Nx, Ny):
    """

    :param c:
    :param Kx:
    :param Ky:
    :param Nx:
    :param Ny:
    :return:
    """
    mesh = MeshGenerator('crazy', c=c, bounds=[(0,1), (0,1)])([Kx, Ky], EDM=None)
    space = SpaceInvoker('polynomials')([('Lobatto', Nx), ('Lobatto', Ny)])
    FC = FormCaller(mesh, space)
    ES = ExactSolutionSelector(mesh)('sL:sincos1')

    u = FC('1-f-o', is_hybrid=False)
    p = FC('2-f-o', is_hybrid=False)
    f = FC('2-f-o', is_hybrid=False)

    B0 = EWC_ColumnVector(mesh, u.num.basis)
    B0.gathering_matrix = u

    f.TW.func.do.set_func_body_as(ES, "source_term")
    f.TW.do.push_all_to_instant(0)
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

    results = LS.solve('GMRES')(0, restart=100, maxiter=20)[0]
    results.do.distributed_to(u, p)

    u.TW.func.do.set_func_body_as(ES, "velocity")
    u.TW.do.push_all_to_instant(0)
    u_error_L2 = u.error.L()
    p.TW.func.do.set_func_body_as(ES, "potential")
    p.TW.do.push_all_to_instant(0)
    p_error_L2 = p.error.L()

    return u_error_L2, p_error_L2