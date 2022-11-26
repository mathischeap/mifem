# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/25 6:19 PM
"""
import sys

import numpy as np

if './' not in sys.path: sys.path.append('/')

from components.freeze.base import FrozenOnly


class mpRfT_SgF_Cochain_Base(FrozenOnly):
    """"""

    def __init__(self, t):
        """"""
        self._t_ = t
        self._trace_ = None
        self._local_ = mpRfT_SgF_Cochain_Base_Local(self)
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
    def trace(self):
        return self._trace_

    @trace.setter
    def trace(self, trace):
        """

        Parameters
        ----------
        trace : dict
            {repr of all segments: the local cochain on segments in a 1-d array.},

        Returns
        -------

        """
        assert isinstance(trace, dict), f"please put local cochain in a dict."
        for rp in self._container_:
            assert rp.__repr__() in trace
        self._trace_ = trace

    @property
    def local(self):
        """Root-Cell-Wise cochain."""
        return self._local_

    def __getitem__(self, rc_rp):
        return self.local[rc_rp]



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
            trace = dict()
            sGM = self._t_.numbering.sgW_gathering
            for seg in sGM:  # go through all local elements
                idx = sGM[seg].full_vector
                trace[seg] = V[idx]
            self.trace = trace

        else:
            raise Exception(f"Can not set cochain from {globe}.")





class mpRfT_SgF_Cochain_Base_Local(FrozenOnly):
    """"""

    def __init__(self, cochain):
        self._cochain_ = cochain
        self._t_ = cochain._t_
        self._mesh_ = self._t_.mesh
        self._freeze_self_()

    def __getitem__(self, rc_rp):
        """Find and return the local cochain in the rc in real time.

        Parameters
        ----------
        rc_rp

        Returns
        -------

        """
        rc = self._mesh_[rc_rp]
        Frame = rc.frame
        local_cochain = list()
        for edge in Frame:
            segments = Frame[edge]
            for seg in segments:
                sgW_cochain = self._cochain_.trace[seg.__repr__()]
                local_cochain.append(sgW_cochain)

        return np.concatenate(local_cochain)








if __name__ == "__main__":
    # mpiexec -n 4 python objects/mpRfT/base/forms/segment/cochain.py
    pass
