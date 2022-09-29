# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 8/13/2022 3:44 PM
"""
import sys
mifem_dir = './' # the dir containing the mifem package
if mifem_dir not in sys.path: sys.path.append(mifem_dir)

from __init__ import tools as mt
from tools.__init__ import linalg
from __init__ import screws as ms
from __init__ import cscg2
from root.config.main import rAnk, mAster_rank
from numpy import pi
import numpy as np
from screws.miscellaneous.mios import rmdir, remove
from screws.miscellaneous.miprint import miprint

def test_Regular_Newton_Raphson():
    """"""
    miprint("RNR [test_Regular_Newton_Raphson] ...... ", flush=True)
    #--------- define the problem ---------------------------------------------------------
    c = 0
    K = 5   # K * K elements (uniform)
    N = 3   # polynomial degree
    dt = 0.01
    t = 0.05
    image_folder = './__images_test_NRR__'
    RDF_filename = 'shear_layer_rollup_p1_NRR_test'
    iterator_name='shear-layer-rollup-p1_NRR_test'
    image_levels = [-6, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 6]
    #--------------------------------------------------------------------------------------
    t0 = 0
    steps = int((t-t0)/dt)
    SI = mt.SimpleIterator(t0=t0, dt=dt, max_steps=steps,
                        auto_save_frequency=10,
                        RDF_filename=RDF_filename,
                        real_time_monitor=True,
                        name=iterator_name)

    IC = ms.Counter()
    ms.mkdir(image_folder)

    mesh = cscg2.mesh('crazy_periodic',
                      bounds=[[0, 2 * pi], [0, 2 * pi]], c=c)(
        [K, K], show_info=True)
    space = cscg2.space('polynomials')([N, N], show_info=False)
    FC = cscg2.form(mesh, space)
    es = cscg2.exact_solution(mesh)("Euler:shear_layer_rollup", show_info=False)

    #----------- unknowns -----------------------------------------------------------------
    u = FC('1-f-o', is_hybrid=False, name='velocity')
    w = FC('0-f-o', is_hybrid=False, name='vorticity')
    P = FC('2-f-o', is_hybrid=False, name='total pressure')

    #--------- tests ----------------------------------------------------------------------
    v = FC('1-f-o', is_hybrid=False, name='test-velocity')
    o = FC('0-f-o', is_hybrid=False, name='test-vorticity')
    q = FC('2-f-o', is_hybrid=False, name='test-total pressure')

    M0 = w.matrices.mass
    M1 = u.matrices.mass
    M2 = P.matrices.mass
    E21 = u.matrices.incidence
    E12 = E21.T
    E10 = w.matrices.incidence
    E01 = E10.T
    C_wk_uk1 = w.special.cross_product_1f__ip_1f(u, v, output='2-M-1')
    C_wk1_uk = w.special.cross_product_1f__ip_1f(u, v, output='2-M-0')
    MDM = w.special.cross_product_1f__ip_1f(u, v, output='MDM')

    Cv = MDM.do.reduce_to_vector(v)

    #----------- initial condition: u, w @ t0 ----------------------------------------------
    u.TW.func.do.set_func_body_as(es, 'velocity')
    u.TW.current_time = t0
    u.TW.do.push_all_to_instant()
    u.discretize()
    KE_t0 = 0.5 * u.do.compute_L2_energy_with(M=M1)
    L2_du_t0 = u.do.compute_Ln_norm_of_coboundary()

    w.TW.func.do.set_func_body_as(es, 'vorticity')
    w.TW.current_time = t0
    w.TW.do.push_all_to_instant()
    w.discretize()
    w.visualize.matplot.contour(levels=image_levels,
                                show_boundaries=False,
                                usetex=False,
                                saveto=image_folder + '/' + str(next(IC)),
                                title=f't=%.3f'%t0)
    EN_t0 = 0.5*w.do.compute_L2_energy_with(M=M0)
    Vor_t0 = w.do.compute_Ln_energy(n=1)

    # ------------ nonlinear system @ tk -------------------------------------------------
    A = ([(1/dt) * M1 + 0.25 * C_wk_uk1, 0.25 * C_wk1_uk, - E12 @ M2],
         [                     E01 @ M1,            - M0,       None],
         [                          E21,            None,       None])

    B = [0.25 * MDM, None, None]

    f = [(1/dt) * M1 @ u - 0.25 * Cv, None, None]

    nLS = linalg.NonLinearSystem([v, o, q], A, B, [u, w, P], f)

    nLS.customize.set_no_evaluation(-1)

    nLS.solve.solver = 'Newton-Raphson'
    nLS.solve.Newton_Raphson.routine = 'regular'

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
        L2_du :
        vorticity :

        """

        assert tk1 == tk + dt, f"A trivial check."
        message = list()

        #---- u, w @ tk -----------------------------------------------------------------------
        x0 = linalg.LocallyFullVector((u, w, P))
        # R = nLS.solve(x0, atol=1e-3, maxiter=5,
                      # LS_solver_para='GMRES', LS_solver_kwargs={'atol':1e-5})
        R = nLS.solve(x0, atol=1e-5, maxiter=5,
                      LS_solver_para='direct', LS_solver_kwargs={})
        R[0].do.distributed_to(u, w, P)
        message.append(R[4])

        # -------- u, w @ tk1 ----------------------------------------------
        w.visualize.matplot.contour(levels=image_levels,
                                    show_boundaries=False,
                                    usetex=False,
                                    saveto=image_folder + '/' + str(next(IC)),
                                    title=f't=%.3f'%tk1)

        EN_tk1  = 0.5 * w.do.compute_L2_energy_with(M=M0)
        KE_tk1  = 0.5 * u.do.compute_L2_energy_with(M=M1)
        Vor_tk1 = w.do.compute_Ln_energy(n=1)
        L2_du_tk1 = u.do.compute_Ln_norm_of_coboundary()

        return 1, 0, message, EN_tk1, KE_tk1, L2_du_tk1, Vor_tk1

    SI(SOLVER, [EN_t0, KE_t0, L2_du_t0, Vor_t0])
    SI.run()

    if rAnk == mAster_rank:
        data = SI.RDF.to_numpy()
        output1 = data[-1][1:]
        output2 = data[-2][1:]
        np.testing.assert_array_almost_equal(output1, output2)

    ms.make_a_video_from_images_in_folder(image_folder, duration=t, clean_images=True)

    remove('shear_layer_rollup_p1_NRR_test.csv')
    remove('MPI_IGR_shear-layer-rollup-p1_NRR_test.png')
    remove(image_folder + '/video.avi')
    rmdir(image_folder)

    return 1



if __name__ == '__main__':
    # mpiexec -n 4 python __tests__/unittests/nonlinear_solver/regular_Newton_Raphson.py
    test_Regular_Newton_Raphson()
