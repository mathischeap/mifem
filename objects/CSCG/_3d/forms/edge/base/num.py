# -*- coding: utf-8 -*-
from components.freeze.base import FrozenOnly



class _3dCSCG_EdgeForm_Num(FrozenOnly):
    """"""
    def __init__(self, ef):
        """"""
        self._ef_ = ef
        self._freeze_self_()


    @property
    def basis(self):
        """(int) Return an int which represent the number of basis function one mesh element has."""
        return self._ef_._NUM_basis_

    @property
    def basis_components(self):
        """Return a dict, keys are corner-edge names, values are the num basis on this corner edge.
        """
        return self._ef_._NUM_basis_components_