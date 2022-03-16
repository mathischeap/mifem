# -*- coding: utf-8 -*-
"""

"""
from screws.freeze.main import FrozenOnly
from _2dCSCG.mesh.domain.boundaries.boundary.main import Boundary


class _2dCSCG_Domain_Boundaries(FrozenOnly):
    """ """
    def __init__(self, domain):
        assert domain.ndim == 2, " <Domain> <Boundaries> "
        assert domain.__class__.__name__ == '_2dCSCG_Domain', " <Domain> <Boundaries> "
        self._domain_ = domain
        self._boundaries_ = {}
        for bn in self.names:
            self._boundaries_[bn] = Boundary(self, bn)
        self._distribution_regularities_ = None
        self._freeze_self_()

    def __getitem__(self, bn):
        """ """
        return self._boundaries_[bn]

    @property
    def ndim(self):
        return 2

    @property
    def num(self):
        return self._domain_._num_boundaries_

    @property
    def names(self):
        return self._domain_._boundary_names_

    @property
    def region_edges(self):
        return self._domain_.domain_input.boundary_region_edges


    @property
    def distribution_regularities(self):
        """How the boundaries are distributed. Return a list containing one or some of:

            (1) to be added...
            (2) ...

        """
        if self._distribution_regularities_ is not None:
            return self._distribution_regularities_

        self._distribution_regularities_ = list()


        #--------- below we do all the regularity checks ---------------------------------

        # TODO: to be added

        #==================================================================================

        return self._distribution_regularities_