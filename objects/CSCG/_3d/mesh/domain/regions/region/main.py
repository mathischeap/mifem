

import numpy as np

from screws.freeze.main import FrozenOnly


from screws.quadrature import Quadrature
from objects.CSCG._3d.mesh.domain.regions.region.interpolations.allocator import InterpolationAllocator
from objects.CSCG._3d.mesh.domain.regions.region.side_geometries.allocator import SideGeometryAllocator
from objects.CSCG._3d.mesh.domain.regions.region.types_wrt_metric.allocator import TypeWr2MetricAllocator
from objects.CSCG._3d.mesh.domain.regions.region.sides.main import Sides

from objects.CSCG._3d.mesh.domain.regions.region.sub_geometry.main import RegionSubGeometry

from objects.CSCG._3d.mesh.domain.regions.region.inheriting.topology import RegionTopology
from objects.CSCG._3d.mesh.domain.regions.region.IS import _3dCSCG_Region_IS


class Region(RegionTopology, FrozenOnly):
    def __init__(self, ndim, name, cc, st, interpolator, domain_input):
        """
        Parameters
        ----------
        ndim : int
            Must be 3.
        name :
        cc : corner coordinates
        st : side type
        interpolator :
        domain_input :
            The DomainInput, but its classname is not 'DomainInput'. It inherits the
            class `DomainInput` but has personal name.

        """
        assert ndim == 3, " < Region> "
        self._ndim_ = ndim
        self._name_ = name
        self._volume_ = None
        self._domain_input_ = domain_input  # can not remove this line
        self.___PRIVATE_set_corner_coordinates___(cc)
        self.___PRIVATE_set_side_types___(st)
        self.___PRIVATE_parse_side_types___()
        self._type_wrt_metric_ = self.___PRIVATE_PARSE_type_wrt_metric___()
        self._interpolation_ = InterpolationAllocator(interpolator)(self)

        self.___is_orthogonal___ = None
        self._sides_ = Sides(self)
        self._sub_geometry_ = RegionSubGeometry(self)
        self._IS_ = None
        self._MAP_ = None
        self._freeze_self_()

    @property
    def map(self):
        return self._MAP_

    @property
    def type_wrt_metric(self):
        return self._type_wrt_metric_

    @property
    def sides(self):
        return self._sides_

    @property
    def name(self):
        return self._name_

    @property
    def ndim(self):
        return self._ndim_

    @property
    def interpolation(self):
        return self._interpolation_

    @property
    def volume(self):
        """ Return the volume of this region."""
        if self._volume_ is None:
            p = 20  # we use this many high order numerical quadrature.
            qx, qy, qz, quad_weights = Quadrature([p, p, p], category='Gauss').quad_ndim_ravel
            detJ = self.interpolation.Jacobian(qx, qy, qz)
            self._volume_ = np.sum(detJ * quad_weights) / 8  # divided by 8 because [-1,1]^3 -> [0,1]^3
        return self._volume_

    @property
    def corner_coordinates(self):
        return self._corner_coordinates_

    @property
    def sub_geometry(self):
        return self._sub_geometry_

    @property
    def IS(self):
        if self._IS_ is None:
            self._IS_ = _3dCSCG_Region_IS(self)
        return self._IS_

    def ___PRIVATE_PARSE_type_wrt_metric___(self):
        if self._domain_input_.region_type_wr2_metric is None:
            return TypeWr2MetricAllocator('chaotic')(self)
        elif self.name in self._domain_input_.region_type_wr2_metric:
            return TypeWr2MetricAllocator(self._domain_input_.region_type_wr2_metric[self.name])(self)
        else:
            return TypeWr2MetricAllocator('chaotic')(self)

    def ___PRIVATE_set_side_types___(self, st):
        """
        Parameters
        ----------
        st : dict
            contain the info of non-plane(3D), non-straight(2D) sides.
        """
        _side_types_ = [None for _ in range(self.num_sides())]
        for key in st:
            if key.split('-')[0] == self.name:
                _side_types_[self._side_name_to_index_(key.split('-')[1])] = st[key]
        self._side_types_ = tuple(_side_types_)

    def ___PRIVATE_parse_side_types___(self):
        """Here we get the 6(3D) side geometries for further getting the region interpolation.

        Attributes
        ----------
        self._side_geometries_ :
            For 3D: N -> S -> W -> E -> B -> F.

        """
        SC_LN = {0: [0, 2, 4, 6],  # the north side has corners locally numbered 0, 2, 4, 6.
                 1: [1, 3, 5, 7],  # the south side has corners locally numbered 1, 3, 5, 7.
                 2: [0, 1, 4, 5],  # the west side has corners locally numbered 0, 1, 4 5.
                 3: [2, 3, 6, 7],  # and so on
                 4: [0, 1, 2, 3],
                 5: [4, 5, 6, 7]}
        _sg_ = dict()
        for i in range(self.num_sides()):
            _st_ = self._side_types_[i]
            if _st_ is None:
                _st_ = ('plane',)
            _cc_ = np.array(self.corner_coordinates)[SC_LN[i]]
            _sg_['NSWEBF'[i]] = SideGeometryAllocator(_st_)(_cc_)
        self._side_geometries_ = _sg_

    def ___PRIVATE_set_corner_coordinates___(self, cc):
        """
        Parameters
        ----------
        cc : tuple
            For 3D:
            np.shape(cc) must be (8, 3), and cc[i] corresponds to NWB, SWB,
            NEB, SEB, NWF, SWF, NEF, SEF corners, and each i is a tuple
            itself represents (x, y, z) coordinates.

        """
        assert np.shape(cc) == (self.num_corners(), self.ndim), \
            " <Region> : coordinates shape={} wrong.".format(np.shape(cc))
        self._corner_coordinates_ = cc