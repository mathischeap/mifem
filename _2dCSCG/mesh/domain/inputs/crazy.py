

from _2dCSCG.mesh.domain.inputs.base import DomainInputBase
import numpy as np


class Crazy(DomainInputBase):
    """"""
    def __init__(self, c=0, bounds=((0, 1), (0, 1))):
        """

        :param c:
        :param bounds:
        """
        assert np.shape(bounds)[0] == 2, " <Crazy> : bounds dimension is wrong."
        for i in range(2):
            assert np.shape(bounds[i]) == (2,) and bounds[i][1] > bounds[i][0], \
                " <CrazyDomain> : bounds[{}]=={} is wrong.".format(i, bounds[i])
        self._bounds_ = bounds
        self._c_ = c
        super().__init__(domain_name='Crazy')
        x0, x1 = bounds[0]
        y0, y1 = bounds[1]
        # _____________ standard inputs ________________________________________________
        self.region_corner_coordinates = {'R:R': ((x0, y0), (x1, y0), (x0, y1), (x1, y1))}
        self.region_edge_types = {'R:R-L': ('free',)}
        self.boundary_region_edges = {'Upper': "R:R-U", 'Down': "R:R-D", 'Left': "R:R-L", 'Right': "R:R-R"}
        self.region_interpolators = 'crazy'

        self.region_type_wr2_metric = {'R:R': 'crazy'}
        self.region_sequence = ('R:R',)
        self.internal_parameters = ['c', ]  # has to be defined after the super().__init__
        # ------------------------------------------------------------------------------

    @property
    def bounds(self):
        return self._bounds_

    @property
    def c(self):
        return self._c_
