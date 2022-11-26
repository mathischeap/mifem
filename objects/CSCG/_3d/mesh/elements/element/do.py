# -*- coding: utf-8 -*-
from components.freeze.base import FrozenOnly
import numpy as np

class _3dCSCG_MeshElement_Do(FrozenOnly):
    """"""
    def __init__(self, element):
        self._element_ = element
        self._freeze_self_()


    def generate_element_plot_data(self, zoom=1, density=10):
        """

        Parameters
        ----------
        zoom : float, default=1
            default The data will be for an element of size `zoom` * [-1,1]^3.
        density : int, default=10

        Returns
        -------
        data :

        """
        mark = self._element_.type_wrt_metric.mark
        if isinstance(mark, str) and mark[:4] == 'Orth':
            density = 2
        else:
            pass

        assert 0 < zoom <= 1, f"zoom={zoom} is wrong!"

        O = np.ones(density) * zoom
        M = - np.ones(density) * zoom
        S = np.linspace(-1 * zoom, 1 * zoom, density)

        data = np.array([
        np.array(self._element_.coordinate_transformation.mapping(S, M, M)),
        np.array(self._element_.coordinate_transformation.mapping(S, M, O)),
        np.array(self._element_.coordinate_transformation.mapping(S, O, O)),
        np.array(self._element_.coordinate_transformation.mapping(S, O, M)),

        np.array(self._element_.coordinate_transformation.mapping(M, S, M)),
        np.array(self._element_.coordinate_transformation.mapping(M, S, O)),
        np.array(self._element_.coordinate_transformation.mapping(O, S, O)),
        np.array(self._element_.coordinate_transformation.mapping(O, S, M)),

        np.array(self._element_.coordinate_transformation.mapping(M, M, S)),
        np.array(self._element_.coordinate_transformation.mapping(M, O, S)),
        np.array(self._element_.coordinate_transformation.mapping(O, O, S)),
        np.array(self._element_.coordinate_transformation.mapping(O, M, S))
        ])

        return data