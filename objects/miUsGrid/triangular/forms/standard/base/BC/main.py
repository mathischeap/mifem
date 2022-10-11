# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 10/4/2022 11:47 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')
from screws.freeze.main import FrozenOnly
from objects.miUsGrid.triangular.forms.standard.base.BC.interpret.main import miUsTriangle_SF_BC_Interpret


class miUsTriangle_SF_BC(FrozenOnly):
    """"""

    def __init__(self, sf):
        """"""
        self._sf_ = sf
        self._CF_ = None
        self._boundaries_ = None
        self._involved_element_edges_ = None
        self._interpret_ = miUsTriangle_SF_BC_Interpret(sf)
        self._freeze_self_()

    @property
    def CF(self):
        return self._CF_

    @CF.setter
    def CF(self, cf):
        self._sf_.___Pr_check_BC_CF___(cf)
        self._CF_ = cf

    @property
    def boundaries(self):
        """The valid boundaries of the BC of this form."""
        return self._boundaries_

    @boundaries.setter
    def boundaries(self, bns):
        """"""

        BNS = self._sf_.mesh.boundaries.names
        if isinstance(bns, str):
            bns = [bns,]
        else:
            pass
        assert isinstance(bns, (list, tuple)), f"pls put boundary names into a list or tuple."
        for i, bn in enumerate(bns):
            assert bn in BNS, f"boundary names [{i}] = {bn} is not a valid boundary name."
        self.___Pr_parse_involved_element_edges___(bns)
        self.interpret.RESET_cache() # reset current interpretations.
        self._boundaries_ = bns


    def ___Pr_parse_involved_element_edges___(self, bns):
        """"""
        RoE = self._sf_.mesh.boundaries.range_of_element_edge
        self._involved_element_edges_ = list()
        for bn in bns:
            self._involved_element_edges_.extend(RoE[bn])

    @property
    def involved_element_edges(self):
        """The involved local element edges."""
        return self._involved_element_edges_

    @property
    def interpret(self):
        return self._interpret_





if __name__ == '__main__':
    # mpiexec -n 4 python objects/miUsGrid/triangular/forms/standard/base/BC/main.py
    from __init__ import miTri
    fc = miTri.form('rand0', 2)
    f0 = fc('0-f-o')

    f0.BC.boundaries = ['Upper', 'Left']

    print(f0.BC.involved_element_edges)