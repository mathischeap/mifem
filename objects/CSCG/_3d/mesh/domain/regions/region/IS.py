# -*- coding: utf-8 -*-

from screws.freeze.base import FrozenOnly




class _3dCSCG_Region_IS(FrozenOnly):
    """"""
    def __init__(self, region):
        """"""
        self._region_ = region
        self._periodic_to_self_ = None
        self._orthogonal_ = None
        self._freeze_self_()

    @property
    def periodic_to_self(self):
        """bool: if this region is periodic to itself.

        The idea is we go through all periodic_boundary_pairs to see if the two boundaries of a
        pair are both in the region map. If yes, then it is periodic to self. Otherwise, it is not.
        """
        if self._periodic_to_self_ is None:
            pbp = self._region_._domain_input_.periodic_boundary_pairs
            MAP = self._region_.map

            T_or_F = False

            for pair in pbp:
                b1, b2 = pair.split('=')
                if b1 in MAP and b2 in MAP:
                    T_or_F = True
                    break

            self._periodic_to_self_ = T_or_F

        return self._periodic_to_self_

    @property
    def orthogonal(self):
        if self._orthogonal_ is None:
            TwM = self._region_.type_wrt_metric
            mark = TwM.mark
            if isinstance(mark, str) and mark[:10] == "orthogonal":
                self._orthogonal_ = True
            elif isinstance(mark, str) and mark[:5] == "crazy" and TwM._c_ == 0:
                    self._orthogonal_ = True
            else:
                self._orthogonal_ = False
        return self._orthogonal_