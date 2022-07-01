# -*- coding: utf-8 -*-


from screws.freeze.main import FrozenOnly


class ____3dCSCG_ADSF_Cochain_Local____(FrozenOnly):
    """"""
    def __init__(self, dsf_CO):
        """"""
        self._PC_ = dsf_CO._dsf_.prime.cochain
        self._MM_ = dsf_CO._dsf_.mass_matrix
        self._freeze_self_()

    def __getitem__(self, i):
        """"""
        return self._MM_[i] @ self._PC_.local[i]

    def __contains__(self, i):
        """"""
        return i in self._PC_.local

    def __iter__(self):
        """Go through all mesh element numbers in this core."""
        for i in self._PC_.local:
            yield i

    def __len__(self):
        """Actually return how many mesh elements in this core."""
        return len(self._PC_.local)