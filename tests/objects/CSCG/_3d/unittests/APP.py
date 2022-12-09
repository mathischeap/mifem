# -*- coding: utf-8 -*-
"""Here we test the programs.

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
import sys
if './' not in sys.path: sys.path.append('./')
from root.config.main import *

from tests.objects.CSCG._3d.unittests.auxiliaries.icpNS_TGV_new import NoHy_TGV_NEW
from tests.objects.CSCG._3d.unittests.auxiliaries.icpsNS_TGV_LGMRES_solver import NoHy_TGV_NEW_LGMRES

import os


def test_APP_NO2_icpsNS_no_hybrid_TGV():
    # warnings.filterwarnings("default")
    if RANK == MASTER_RANK:
        print("--- [test_APP_NO2_icpsNS_no_hybrid_TGV] ...... ", flush=True)

    N, k, t, steps, Re = 2, 3, 1, 5, 500


    SN = NoHy_TGV_NEW(N=N, k=k, t=t, steps=steps, Re=Re, atol=1e-5, restart=50, maxiter=10, show_info=False, save_uw=False)

    SL = NoHy_TGV_NEW_LGMRES(N=N, k=k, t=t, steps=steps, Re=Re, atol=1e-5, m=45, K=5, maxiter=10, show_info=False, save_uw=False)

    if RANK == MASTER_RANK:
        os.remove(f'icpNS_NH_RDF___TGV_NEW_Re{Re}_N{N}K{k}t{t}Steps{steps}.csv')
        os.remove(f'icpNS_NH_RDF_LGMRES__TGV_NEW_Re{Re}_N{N}K{k}t{t}Steps{steps}.csv')

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
    # mpiexec -n 4 python tests/objects/CSCG/_3d/unittests/APP.py
    test_APP_NO2_icpsNS_no_hybrid_TGV()