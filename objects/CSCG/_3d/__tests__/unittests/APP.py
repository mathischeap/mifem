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
from root.config.main import *

from objects.CSCG._3d.__tests__.unittests.auxiliaries.icpsNS_manu_conserving_solver import manu_conserving_solver
from objects.CSCG._3d.__tests__.unittests.auxiliaries.icpsNS_TGV_no_hybrid import NoHy_TGV
from objects.CSCG._3d.__tests__.unittests.auxiliaries.icpNS_TGV_new import NoHy_TGV_NEW
import os
import random
import warnings
from objects.CSCG._3d.__tests__.unittests.auxiliaries.icpsNS_TGV_LGMRES_solver import NoHy_TGV_NEW_LGMRES





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

    N, k, t, steps, Re = 2, 3, 1, 5, 500

    SI = NoHy_TGV(N=N, k=k, t=t, steps=steps, Re=Re, tol=1e-5, restart=50, maxiter=10, show_info=False, save_uw=False)

    SN = NoHy_TGV_NEW(N=N, k=k, t=t, steps=steps, Re=Re, atol=1e-5, restart=50, maxiter=10, show_info=False, save_uw=False)

    SL = NoHy_TGV_NEW_LGMRES(N=N, k=k, t=t, steps=steps, Re=Re, atol=1e-5, m=45, K=5, maxiter=10, show_info=False, save_uw=False)

    if rAnk == mAster_rank:
        os.remove(f'RDF_TGV_Re{Re}_N{N}K{k}t{t}Steps{steps}.csv')
        os.remove(f'icpNS_NH_RDF___TGV_NEW_Re{Re}_N{N}K{k}t{t}Steps{steps}.csv')
        os.remove(f'icpNS_NH_RDF_LGMRES__TGV_NEW_Re{Re}_N{N}K{k}t{t}Steps{steps}.csv')

        data = SI.RDF.to_numpy()
        np.testing.assert_array_almost_equal(data[-1], np.array(
            [1.00000000e+00, 2.00000000e-01, 1.13837441e-01, 1.13745439e-01, 1.16437118e-01,
             -1.23027780e-05, -1.48929594e-05, 3.93532257e-01, 4.05023500e-01, 4.17056613e-01,
             2.55458022e+00, 6.47060331e+00, 3.13155199e-04]), decimal=5)

        data = SN.RDF.to_numpy()
        np.testing.assert_array_almost_equal(data[-1], np.array(
            [1.00000000e+00, 2.00000000e-01, 1.13837352e-01, 1.13745175e-01, 1.16437328e-01, 7.36504537e-08,
             7.08131196e-08, 3.93534048e-01, 4.05025338e-01, 4.17058241e-01, 2.55462821e+00, 6.47077101e+00,
             1.19959342e-06]), decimal=5)

        data = SL.RDF.to_numpy()
        np.testing.assert_array_almost_equal(data[-1], np.array(
            [1.00000000e+00, 2.00000000e-01, 1.13837352e-01, 1.13745175e-01, 1.16437328e-01,
             -2.76292038e-08, -2.77394652e-08, 3.93534038e-01, 4.05025339e-01, 4.17058242e-01,
             2.55462822e+00, 6.47077097e+00, 2.29692361e-09]), decimal=5)

    return 1






if __name__ == '__main__':
    # mpiexec -n 4 python objects\CSCG\_3d\__tests__\unittests\APP.py
    test_APP_NO1_icpsNS_no_hybrid_manu_conserving()
    test_APP_NO2_icpsNS_no_hybrid_TGV()