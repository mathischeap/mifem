# -*- coding: utf-8 -*-
from components.freeze.base import FrozenOnly



class ElementPerpendicularSlice(FrozenOnly):
    """An element slice can only be a full perpendicular surface (refer to the reference element,
    can be curvilinear in the physical element.)

    For example, perpendicular_to_axis = 'xi', position = 0, then it is the surface xi=0, eta=sigma=[-1,1].
    """

    def __init__(self, element, xi=None, eta=None, sigma=None):
        self._element_ = element

        NUM_NONE = 0
        for j, _ in enumerate([xi, eta, sigma]):
            if _ is None:
                NUM_NONE += 1
            else:
                perpendicular_to_axis = j
                assert isinstance(_, (int, float)) and -1 <= _ <= 1, f"position = {_} is wrong."
                self._position_ = _

        assert NUM_NONE == 2, f'xi={xi}, eta= {eta}, sigma= {sigma} wrong.'
        # noinspection PyUnboundLocalVariable
        self._PTA_ = ['xi', 'eta', 'sigma'][perpendicular_to_axis]

        self._freeze_self_()


    @property
    def perpendicular_to_axis(self):
        """"""
        return self._PTA_


    @property
    def position(self):
        """A value between -1 and 1 to indicate the slice location."""
        return self._position_