# -*- coding: utf-8 -*-


from components.freeze.base import FrozenOnly



class _3dCSCG_1LocalTrace_Discretize_StandardVector_T_perp(FrozenOnly):
    """"""
    def __init__(self, tf):
        self._tf_ = tf
        self._freeze_self_()

    def __call__(self,
        update_cochain=True, target='func'):
        """We will discretize the Trace_perpendicular component of a standard vector field to all trace
        elements.

        :param update_cochain:
        :param target:
        :return:
        """
        if target in ('BC',): assert update_cochain is False, f"CANNOT update cochain when target is {target}"
        raise NotImplementedError()