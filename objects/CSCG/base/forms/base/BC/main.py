# -*- coding: utf-8 -*-
from components.freeze.base import FrozenOnly
from objects.CSCG.base.forms.base.BC.interpret.main import CSCG_FORM_BC_Interpret


class CSCG_Form_BC(FrozenOnly):
    def __init__(self, f):
        self._f_ = f
        self._CF_ = None
        self._boundaries_ = list()
        self._involved_element_parts_ = None
        self._EWC_ = None
        self._interpret_ = CSCG_FORM_BC_Interpret(self._f_)
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

        if bns == 'all':
            bns = BNS
        else:
            pass

        if isinstance(bns, str):
            bns = [bns, ]
        else:
            pass
        assert isinstance(bns, (list, tuple)), f"pls put boundary names into a list or tuple."
        for i, bn in enumerate(bns):
            assert bn in BNS, f"boundary names [{i}] = {bn} is not a valid boundary name."
        self.___Pr_parse_involved_element_parts___(bns)
        self._boundaries_ = bns

    def ___Pr_parse_involved_element_parts___(self, bns):
        """"""
        mesh = self._f_.mesh
        if mesh.ndim == 3:
            Res = mesh.boundaries.range_of_element_sides
        elif mesh.ndim == 2:
            Res = mesh.boundaries.range_of_element_edges
        else:
            raise Exception()

        self._involved_element_parts_ = list()
        for bn in bns:
            self._involved_element_parts_.extend(Res[bn])

    @property
    def interpret(self):
        """Use no cache in interpret, it will generate data in real time.
        Remember to save its sub-properties with separate variables.
        """
        return self._interpret_
