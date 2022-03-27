from screws.freeze.base import FrozenOnly


class _2dCSCG_Domain_IS(FrozenOnly):
    """"""

    def __init__(self, domain):
        """"""
        self._d_ = domain
        self._freeze_self_()

    @property
    def periodic(self):
        """This domain involves periodic boundary."""
        return len(self._d_.domain_input.periodic_boundary_pairs) != 0

    @property
    def fully_periodic(self):
        """This domain is fully periodic, so \partial \Omega = \empty"""
        return self._d_.domain_input.periodic_boundaries == set(self._d_.boundaries.names)