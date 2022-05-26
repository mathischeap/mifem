# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/19 6:29 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

import numpy as np
from screws.freeze.base import FrozenOnly


class mpRfT_SF_Cochain_Base(FrozenOnly):
    """"""

    def __init__(self, f):
        """"""
        self._f_ = f
        self._local_ = None
        self.___Pr_reset_cache___()
        self._freeze_self_()

    def ___Pr_reset_cache___(self): # when we have new representations of cochain, add their caches here
        """"""

    def ___Pr_RGW_cochain___(self):
        raise NotImplementedError()


    @property
    def local(self):
        return self._local_

    @local.setter
    def local(self, local):
        """

        Parameters
        ----------
        local : dict
            {repr of all local-root-cells: the local cochain in a 1-d array.},

        Returns
        -------

        """
        assert isinstance(local, dict), f"please put local cochain in a dict."
        for rp in self._f_.mesh.rcfc:
            lcc = local[rp]
            num_basis = self._f_.num.basis[rp]
            assert lcc.shape == (num_basis,), f"lcc shape dis-match for cell:{rp}"
        self._local_ = local


    def ___Pr_divide_local___(self, rp):
        """"""
        local_rp = self.local[rp]
        return np.split(local_rp, self._f_.num.components)






if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
