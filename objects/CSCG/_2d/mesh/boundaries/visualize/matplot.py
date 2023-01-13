# -*- coding: utf-8 -*-
from components.freeze.base import FrozenOnly
import numpy as np


class _2dCSCG_Mesh_Boundaries_Matplot(FrozenOnly):
    """"""
    def __init__(self, boundaries):
        """"""
        self._boundaries_ = boundaries
        self._freeze_self_()

    def __call__(self, *args, **kwargs):
        return self.___outline___(*args, **kwargs)

    def ___outline___(
            self,
            density=20,
            data_only=False,
    ):
        """

        Parameters
        ----------
        data_only

        Returns
        -------

        """

        mesh = self._boundaries_._mesh_
        RegionEdges = mesh.domain.boundaries.region_edges
        Regions = mesh.domain.regions

        ONES = np.ones(density)
        coo = {
            'U': [0 * ONES, np.linspace(0, 1, density)],
            'D': [ONES, np.linspace(0, 1, density)],
            'L': [np.linspace(0, 1, density), 0 * ONES],
            'R': [np.linspace(0, 1, density), ONES],
        }

        outlines = dict()

        for bn in RegionEdges:
            outlines[bn] = list()
            region_edges = RegionEdges[bn]
            for region_edge in region_edges:
                region, edge = region_edge.split('-')

                region = Regions[region]

                interpolation = region.interpolation

                otl = interpolation.mapping(*coo[edge])

                outlines[bn].append(otl)

        if data_only:
            return outlines
        else:
            raise NotImplementedError()
