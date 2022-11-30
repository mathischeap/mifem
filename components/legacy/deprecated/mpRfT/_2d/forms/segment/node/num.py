# -*- coding: utf-8 -*-
"""
INTRO

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
import sys

if './' not in sys.path: sys.path.append('/')

from components.freeze.base import FrozenOnly


class mpRfT2_NSgF_Num(FrozenOnly):
    """"""

    def __init__(self, t):
        """"""
        self._t_ = t
        self._basis_ = ___Pr_Basis___(t)
        self._freeze_self_()

    @property
    def basis(self):
        return self._basis_

    @property
    def local_dofs(self):
        if self._t_.numbering._num_local_dofs_ is None:
            _ = self._t_.numbering.sgW_gathering
        return self._t_.numbering._num_local_dofs_

    @property
    def GLOBAL_dofs(self):
        return self._t_.numbering.gathering.global_num_dofs


class ___Pr_Basis___(FrozenOnly):
    """"""
    def __init__(self, t):
        self._t_ = t
        self._mesh_ = self._t_.mesh
        self._freeze_self_()

    def __getitem__(self, seg_or_rc_rp):
        """"""
        if seg_or_rc_rp.__class__.__name__ == 'mpRfT2_Segment':
            N = self._t_.N[seg_or_rc_rp]
            return N + 1
        elif isinstance(seg_or_rc_rp, str):
            cell = self._mesh_[seg_or_rc_rp]
            frame = cell.frame
            N = 0
            for edge in frame:
                segments = frame[edge]
                for seg in segments:
                    N += self._t_.N[seg] + 1
            return N




