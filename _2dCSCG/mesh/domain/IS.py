from screws.freeze.inheriting.frozen_only import FrozenOnly


class _2dCSCG_Domain_IS(FrozenOnly):
    """"""

    def __init__(self, domain):
        """"""
        self._d_ = domain
        self._freeze_self_()

    @property
    def periodic(self):
        return True if len(self._d_.domain_input.periodic_boundary_pairs) != 0 else False

