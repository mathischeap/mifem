# -*- coding: utf-8 -*-
from components.freeze.base import FrozenOnly




class _3dCSCG_Domain_IS(FrozenOnly):
    """"""
    def __init__(self, domain):
        """"""
        self._d_ = domain
        self._periodic_ = None
        self._fully_periodic_ = None
        self._freeze_self_()

    @property
    def periodic(self):
        """If this computation domain involve periodic boundaries."""
        if self._periodic_ is None:
            self._periodic_ = len(self._d_.domain_input.periodic_boundary_pairs) != 0
        return self._periodic_

    @property
    def fully_periodic(self):
        """This domain is fully periodic, so \partial \Omega = \empty"""
        if self._fully_periodic_ is None:
            self._fully_periodic_ = self._d_.domain_input.periodic_boundaries == \
                                    set(self._d_.boundaries.names)
        return self._fully_periodic_