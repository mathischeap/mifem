# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/25 6:19 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly


class mpRfT_SgF_Cochain_Base(FrozenOnly):
    """"""

    def __init__(self, t):
        """"""
        self._t_ = t
        self._local_ = None
        self.___Pr_reset_cache___()
        if t.ndim == 2:
            self._container_ = t.mesh.segments
        elif t.ndim == 3:
            self._container_ = t.mesh.facets
        else:
            raise Exception()

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
            {repr of all segments: the local cochain on segments in a 1-d array.},

        Returns
        -------

        """
        assert isinstance(local, dict), f"please put local cochain in a dict."
        for rp in self._container_:
            assert rp.__repr__() in local
        self._local_ = local


if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
