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

from components.miscellaneous.mios import remove, isfile

def test_APP_NO2_icpsNS_no_hybrid_TGV():
    # warnings.filterwarnings("default")
    if RANK == MASTER_RANK:
        print("--- [test_APP_NO2_icpsNS_no_hybrid_TGV] ...... ", flush=True)

    N, k, t, steps, Re = 2, 3, 1, 5, 500

    filename1 = f'icpNS_NH_RDF___TGV_NEW_Re{Re}_N{N}K{k}t{t}Steps{steps}.csv'
    filename2 = f'icpNS_NH_RDF_LGMRES__TGV_NEW_Re{Re}_N{N}K{k}t{t}Steps{steps}.csv'

    if isfile(filename1): remove(filename1)
    if isfile(filename2): remove(filename2)

    SN = NoHy_TGV_NEW(N=N, k=k, t=t, steps=steps, Re=Re, atol=1e-5, restart=50, maxiter=10, show_info=False, save_uw=False)

    SL = NoHy_TGV_NEW_LGMRES(N=N, k=k, t=t, steps=steps, Re=Re, atol=1e-5, m=45, K=5, maxiter=10, show_info=False, save_uw=False)

    if RANK == MASTER_RANK:
        os.remove(filename1)
        os.remove(filename2)

        data = SN.RDF.to_numpy()
        np.testing.assert_array_almost_equal(data[-1], np.array(
            [1.0000e+00, 2.0000e-01, 1.1384e-01, 1.1375e-01, 1.1644e-01,
             7.3639e-08, 7.0782e-08, 3.9353e-01, 4.0503e-01, 4.1706e-01,
             2.5546e+00, 6.4708e+00, 1.1996e-06]), decimal=4)

        data = SL.RDF.to_numpy()
        np.testing.assert_array_almost_equal(data[-1], np.array(
            [1.0000e+00, 2.0000e-01, 1.1384e-01, 1.1375e-01, 1.1644e-01,
             -2.7637e-08, -2.7747e-08, 3.9353e-01, 4.0503e-01, 4.1706e-01,
             2.5546e+00, 6.4708e+00, 2.2965e-09]), decimal=4)

    return 1






if __name__ == '__main__':
    # mpiexec -n 4 python tests/objects/CSCG/_3d/unittests/APP.py
    test_APP_NO2_icpsNS_no_hybrid_TGV()