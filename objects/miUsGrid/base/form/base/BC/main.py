# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 10/4/2022 11:47 PM
"""
from components.freeze.main import FrozenOnly
from objects.miUsGrid.base.form.base.BC.interpret.main import miUsForm_Form_BC_Interpret


class miUsGrid_Form_BC(FrozenOnly):
    """"""

    def __init__(self, f):
        """"""
        self._f_ = f
        self._CF_ = None
        self._boundaries_ = None
        self._involved_element_parts_ = None
        self._interpret_ = miUsForm_Form_BC_Interpret(self._f_)
        self._freeze_self_()

    @property
    def CF(self):
        return self._CF_

    @CF.setter
    def CF(self, cf):
        if cf is None:
            self._CF_ = None
        else:
            self._f_.___Pr_check_BC_CF___(cf)
            self._CF_ = cf

    @property
    def boundaries(self):
        """The valid boundaries of the BC of this form."""
        return self._boundaries_

    @boundaries.setter
    def boundaries(self, bns):
        """"""

        BNS = self._f_.mesh.boundaries.names
        if isinstance(bns, str):
            bns = [bns,]
        else:
            pass
        assert isinstance(bns, (list, tuple)), f"pls put boundary names into a list or tuple."
        for i, bn in enumerate(bns):
            assert bn in BNS, f"boundary names [{i}] = {bn} is not a valid boundary name."
        self.___Pr_parse_involved_element_parts___(bns)
        self._boundaries_ = bns

    def ___Pr_parse_involved_element_parts___(self, bns):
        """"""
        if self._f_.ndim == 2:
            RoE = self._f_.mesh.boundaries.range_of_element_edge
        else:
            raise NotImplementedError()

        self._involved_element_parts_ = list()
        for bn in bns:
            self._involved_element_parts_.extend(RoE[bn])

    @property
    def interpret(self):
        """Use no cache in interpret, it will generate data in real time.
        Remember to save its sub-properties with separate variables.
        """
        return self._interpret_