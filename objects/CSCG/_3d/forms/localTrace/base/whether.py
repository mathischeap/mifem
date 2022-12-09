# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 12/2/2022 9:12 PM
"""
from components.freeze.main import FrozenOnly

class _3dCSCG_LocalTrace_Whether(FrozenOnly):
    """"""

    def __init__(self, ltf):
        """"""
        self._ltf_ = ltf
        self._freeze_self_()

    @property
    def hybrid(self):
        """If it is hybrid in mesh-elements locally?

        Note that, across mesh-elements, local-trace forms must be hybrid. We have no property
        referring to that.
        """
        return self._ltf_._hybrid_

    @property
    def inner_oriented(self):
        return True if self._ltf_.orientation == 'inner' else False