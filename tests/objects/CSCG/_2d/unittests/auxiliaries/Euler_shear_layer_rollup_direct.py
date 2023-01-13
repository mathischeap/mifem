# -*- coding: utf-8 -*-
"""the second program for the 2d inviscid shear lay rollup test case, using direct solver"""

from numpy import pi
from root.config.main import RANK, MASTER_RANK
import os

from tools.iterators.simple import SimpleIterator
from objects.CSCG._2d.master import MeshGenerator, SpaceInvoker, FormCaller, ExactSolutionSelector
from tools.elementwiseCache.dataStructures.objects.columnVector.main import EWC_ColumnVector
from tools.elementwiseCache.dataStructures.operators.bmat.main import bmat
from tools.elementwiseCache.dataStructures.operators.concatenate.main import concatenate
from tools.miLinearAlgebra.linearSystem.main import LinearSystem
from components.generators.counter import Counter

from components.video.make_a_video_from_images_in_folder import make_a_video_from_images_in_folder


def Euler_shear_layer_rollup_direct_test(K, N, dt, t, image_folder, RDF_filename):
    # K = 20 # K * K elements (uniform)
    # N = 2  # polynomial degree
    # dt = 0.01
    # t = 8
    # image_folder = './images_direct'
    # RDF_filename = 'shear_layer_rollup_p2_direct'
    # iterator_name='shear-layer-rollup-p2-direct'
    image_levels = [-6, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 6]

    t0 = 0
    steps = int((t-t0)/dt)
    SI = SimpleIterator(t0=t0, dt=dt, max_steps=steps,
                        RDF_filename=RDF_filename,
                        name='iterator_name')
    IC = Counter()
    if RANK == MASTER_RANK:
        if os.path.isdir(image_folder):
            pass
        else:
            os.mkdir(image_folder)

    mesh = MeshGenerator('rectangle_periodic',
                         p_UL=(0, 0), width=2 * pi, length=2 * pi,
                         region_layout=(2, 2))([K, K])

    space = SpaceInvoker('polynomials')([('Lobatto', N), ('Lobatto', N)])
    FC = FormCaller(mesh, space)
    es = ExactSolutionSelector(mesh)("Euler:shear_layer_rollup")

    w = FC('0-f-o', hybrid=False, name='vorticity')
    u = FC('1-f-o', hybrid=False, name='velocity')
    P = FC('2-f-o', hybrid=False, name='total pressure')

    M0 = w.matrices.mass
    M1 = u.matrices.mass
    M2 = P.matrices.mass
    E21 = u.matrices.incidence
    E12 = E21.T
    C = w.special.cross_product_1f__ip_1f(u, u)

    # ---------- initial condition -----------------------------------------------------------
    u.CF = es.velocity
    u.CF.current_time = t0
    u.discretize()
    KE_t0 = 0.5 * u.do.compute_L2_energy_with(M=M1)
    L2_du_t0 = u.do.compute_Ln_norm_of_coboundary()

    w.CF = es.vorticity
    w.CF.current_time = t0
    w.discretize()
    w.visualize.matplot.contour(levels=image_levels,
                                show_boundaries=False,
                                saveto=image_folder + '/' + str(next(IC)),
                                title=f't=%.3f' % t0)

    # --------- 1/2 step ----------------------------------------------------------------

    A = ([(2/dt) * M1, - E12 @ M2],
         [E21,   None])
    A = bmat(A)
    A.gathering_matrices = ((u, P), (u, P))

    b1 = ((2/dt) * M1 - C) @ u
    b2 = EWC_ColumnVector(mesh, P.num.basis)
    b = concatenate([b1, b2])
    b.gathering_matrix = (u, P)

    LS = LinearSystem(A, b)
    LS.customize.identify_global_row(-1)

    result = LS.solve('direct')()
    u_t0_half = FC('1-f-o', hybrid=False, name='velocity_t0_half')
    result[0].do.distributed_to(u_t0_half, P)

    E01 = w.matrices.incidence.T
    b = concatenate([E01 @ M1 @ u_t0_half, ])
    A = bmat(([M0, ],))
    A.gathering_matrices = (w, w)
    b.gathering_matrix = w
    LS = LinearSystem(A, b)
    result = LS.solve('direct')()
    result[0].do.distributed_to(w)

    EN_t0_half = 0.5 * w.do.compute_L2_energy_with(M=M0)
    Vor_t0_half = w.do.compute_Ln_energy(n=1)

    # ------------ u @ t0,   w @ t_1/2 --------------------------------------------------
    A = ([(1/dt) * M1 + 0.5 * C, - E12 @ M2],
         [E21,   None])
    A = bmat(A)
    A.gathering_matrices = ((u, P), (u, P))

    b1 = ((1/dt) * M1 - 0.5 * C) @ u.cochain.EWC
    b2 = EWC_ColumnVector(mesh, P.num.basis)
    b = concatenate([b1, b2])
    b.gathering_matrix = (u, P)

    LSuP = LinearSystem(A, b)
    LSuP.customize.identify_global_row(-1)
    LSuP.A.do.lock_sparsity()

    A = bmat(([0.5 * M0, ],))
    A.gathering_matrices = (w, w)

    b = concatenate([E01 @ M1 @ u - 0.5 * M0 @ w, ])
    b.gathering_matrix = w
    LS_w = LinearSystem(A, b)
    LS_w.A.do.lock_sparsity()
    LS_w.A.do.lock_assembled_matrix()

    w.visualize.matplot.contour(levels=image_levels,
                                show_boundaries=False,
                                saveto=image_folder + '/' + str(next(IC)),
                                title=f't=%.3f' % (dt/2))

    def SOLVER(tk, tk1):
        """
        Parameters
        ----------
        tk :
        tk1 :

        Returns
        -------
        exit_code: The standard exit code.
        shut_down: If it is ``True``, the outer iterator will shut down immediately.
        message: The solver message.
        enstrophy :
        kinetic_energy :
        vorticity :
        L2_du:

        """

        assert tk1 == tk + dt, f"A trivial check."
        message = list()

        # --- u @ tk , w @ tk+half
        R = LSuP.solve('direct')()
        R[0].do.distributed_to(u, P)
        message.append(R[4])
        # ------- u @ tk+1 ---------------------
        R = LS_w.solve('direct')()
        R[0].do.distributed_to(w)
        message.append(R[4])
        # -------- w @ tk1+half --------------------------------
        w.visualize.matplot.contour(levels=image_levels,
                                    show_boundaries=False,
                                    saveto=image_folder + '/' + str(next(IC)),
                                    title=f't=%.3f' % (tk1+dt/2))
        EN_tk1_half = 0.5 * w.do.compute_L2_energy_with(M=M0)
        KE_tk1 = 0.5 * u.do.compute_L2_energy_with(M=M1)
        Vor_tk1_half = w.do.compute_Ln_energy(n=1)
        L2_du_tk1 = u.do.compute_Ln_norm_of_coboundary()

        return 1, 0, message, EN_tk1_half, KE_tk1, Vor_tk1_half, L2_du_tk1

    SI(SOLVER, [EN_t0_half, KE_t0, Vor_t0_half, L2_du_t0])
    SI.run()

    make_a_video_from_images_in_folder(image_folder, duration=t, clean_images=True)

    return SI
