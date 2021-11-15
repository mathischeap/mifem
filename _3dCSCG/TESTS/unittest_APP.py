# -*- coding: utf-8 -*-
"""
Here we test the programs.

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
import sys
if './' not in sys.path: sys.path.append('./')
from root.config import *
from _3dCSCG.APP.contents.icpsNS.no_hybrid.manu_conserving import manu_conserving_solver
from _3dCSCG.APP.contents.icpsNS.no_hybrid.TGV import NoHy_TGV
from _3dCSCG.APP.contents.icpsNS.no_hybrid.TGV_new import NoHy_TGV_NEW
import os
import random
import warnings
from _3dCSCG.TESTS.__unittest_scripts__.icpsNS_TGV_LGMRES_solver import NoHy_TGV_NEW_LGMRES

def test_APP_NO1_icpsNS_no_hybrid_manu_conserving():
    if rAnk == mAster_rank:
        print(">>> [test_APP_NO1_icpsNS_no_hybrid_manu_conserving] ...... ", flush=True)

    if rAnk == mAster_rank:
        i = random.randint(0, 1)
        s = random.randint(5, 6)
        if i == 0:
            N, k, t, steps = 2, 3, 1, s
        elif i == 1:
            N, k, t, steps = 3, 2, 1, s
        else:
            raise Exception()
    else:
        N, k, t, steps = None, None, None, None
    N, k, t, steps = cOmm.bcast([N, k, t, steps], root=mAster_rank)

    warnings.filterwarnings("ignore")
    SI = manu_conserving_solver(N, k, t, steps)
    warnings.filterwarnings("default")

    if rAnk == mAster_rank:
        os.remove(f'RDF_manuC_N{N}K{k}t{t}Steps{steps}.csv')

        data = SI.RDF.to_numpy()
        output0 = data[1]
        output1 = data[-1]
        output = output1 - output0

        np.testing.assert_almost_equal(output[3], 0)
        np.testing.assert_almost_equal(output[4], 0)
        np.testing.assert_almost_equal(output[5], 0)
        np.testing.assert_almost_equal(output[6], 0)


    return 1


def test_APP_NO2_icpsNS_no_hybrid_TGV():

    # warnings.filterwarnings("default")

    if rAnk == mAster_rank:
        print("--- [test_APP_NO2_icpsNS_no_hybrid_TGV] ...... ", flush=True)

    N, k, t, steps, Re = 2, 2, 1, 5, 500

    SI = NoHy_TGV(N=N, k=k, t=t, steps=steps, Re=Re, tol=1e-5, restart=50, maxiter=10, show_info=False, save_uw=False)

    SN = NoHy_TGV_NEW(N=N, k=k, t=t, steps=steps, Re=Re, atol=1e-5, restart=50, maxiter=10, show_info=False, save_uw=False)

    SL = NoHy_TGV_NEW_LGMRES(N=N, k=k, t=t, steps=steps, Re=Re, atol=1e-5, m=45, K=5, maxiter=10, show_info=False, save_uw=False)

    if rAnk == mAster_rank:
        os.remove(f'RDF_TGV_Re{Re}_N{N}K{k}t{t}Steps{steps}.csv')
        os.remove(f'icpNS_NH_RDF___TGV_NEW_Re{Re}_N{N}K{k}t{t}Steps{steps}.csv')
        os.remove(f'icpNS_NH_RDF_LGMRES__TGV_NEW_Re{Re}_N{N}K{k}t{t}Steps{steps}.csv')

        data = SI.RDF.to_numpy()
        np.testing.assert_array_almost_equal(data[-1], np.array(
            [1.00000000e+00, 2.00000000e-01, 4.43469351e-02, 4.42921795e-02,
             1.53859088e-01, 4.59871030e-11, 1.12474918e-10, 4.67675548e-01,
             1.74664933e-01, 1.77093370e-01, 4.98026933e+00, 9.75551514e+00,
             4.26487127e-06]), decimal=5)

        data = SN.RDF.to_numpy()
        np.testing.assert_array_almost_equal(data[-1], np.array(
            [1.00000000e+00, 2.00000000e-01, 4.43468559e-02, 4.42920630e-02,
             1.53859088e-01, 4.78745965e-09, 4.78746109e-09, 4.67675548e-01,
             1.74664643e-01, 1.77092951e-01, 4.98027219e+00, 9.75551833e+00,
             2.93994455e-15]), decimal=5)

        data = SL.RDF.to_numpy()
        np.testing.assert_array_almost_equal(data[-1], np.array(
            [1.00000000e+00, 2.00000000e-01, 4.43468559e-02, 4.42920630e-02,
             1.53859088e-01, 4.78745965e-09, 4.78746109e-09, 4.67675548e-01,
             1.74664643e-01, 1.77092951e-01, 4.98027219e+00, 9.75551833e+00,
             2.93994455e-15]), decimal=5)


    return 1


if __name__ == '__main__':
    # mpiexec python _3dCSCG\TESTS\unittest_APP.py
    test_APP_NO1_icpsNS_no_hybrid_manu_conserving()
    test_APP_NO2_icpsNS_no_hybrid_TGV()