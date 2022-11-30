# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/19 6:29 PM
"""
import sys

if './' not in sys.path: sys.path.append('/')

import numpy as np
from scipy.sparse import csr_matrix
from components.freeze.base import FrozenOnly
from tools.elementwiseCache.dataStructures.objects.columnVector.main import EWC_ColumnVector


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

    def __getitem__(self, rc_rp):
        return self.local[rc_rp]

    def ___Pr_divide_local___(self, rp):
        """"""
        local_rp = self.local[rp]
        return np.split(local_rp, self._f_.num.components)

    @property
    def EWC(self):
        """Return the cochain as an EWC_ColumnVector."""
        ewc = EWC_ColumnVector(self._f_.mesh, self.___Pr_EWC_DG___, cache_key_generator = 'no_cache')
        ewc.gathering_matrix = self._f_
        return ewc

    def ___Pr_EWC_DG___(self, rc_rp):
        """"""
        return csr_matrix(self.local[rc_rp]).T


    @property
    def globe(self):
        raise NotImplementedError()

    @globe.setter
    def globe(self, globe):
        """
        This process is complex, but it makes sure that the distribution is correct for all cases.

        :param globe:
        :return:
        """
        if globe.__class__.__name__ == 'LocallyFullVector':
            V = globe.V # V already be 1-d array.
            local = dict()
            GM = self._f_.numbering.gathering
            for rc_rp in GM:  # go through all local elements
                idx = GM[rc_rp].full_vector
                local[rc_rp] = V[idx]
            self.local = local

        else:
            raise Exception(f"Can not set cochain from {globe}.")

if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
