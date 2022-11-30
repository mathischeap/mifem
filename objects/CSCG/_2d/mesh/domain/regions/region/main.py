# -*- coding: utf-8 -*-


from objects.CSCG._2d.mesh.domain.regions.region.interpolations.allocator import InterpolationSearcher
from objects.CSCG._2d.mesh.domain.regions.region.edge_geometries.allocator import EdgeGeometryDispatcher
from objects.CSCG._2d.mesh.domain.regions.region.types_wrt_metric.allocator import TypeWr2MetricGiver
from objects.CSCG._2d.mesh.domain.regions.region.whether import _2dCSCG_Region_Whether

from components.freeze.main import FrozenOnly

from objects.CSCG._2d.mesh.domain.regions.region.edges.main import Edges

from components.decorators.all import accepts

from objects.CSCG._2d.mesh.domain.regions.region.inheriting.topology import RegionTopology
import numpy as np




class Region(RegionTopology, FrozenOnly):
    @accepts('self', int, str)
    def __init__(self, ndim, name, cc, et, interpolator, domain_input):
        """
        Parameters
        ---------
        ndim : int
            Must be 2.
        name :
        cc : corner coordinates
        et :
        interpolator :
        domain_input :
            The DomainInput, but its classname is not 'DomainInput'. It inherits the
            class `DomainInput` but has personal name.

        """
        assert ndim == 2, " < Region> "
        self._ndim_ = ndim
        self._name_ = name
        self._domain_input_ = domain_input  # can not remove this line
        self.___PRIVATE_set_corner_coordinates___(cc)
        self.___PRIVATE_set_edge_types___(et)
        self.___PRIVATE_parse_edge_types___()
        self._type_wrt_metric_ = self.___PRIVATE_PARSE_type_wrt_metric___()
        self._interpolation_ = InterpolationSearcher(interpolator)(self)
        self._edges_ = Edges(self)
        self._MAP_ = None
        self._whether_ = None
        self._freeze_self_()

    @property
    def map(self):
        return self._MAP_

    @property
    def whether(self):
        if self._whether_ is None:
            self._whether_ = _2dCSCG_Region_Whether(self)
        return self._whether_

    @property
    def edges(self):
        return self._edges_

    @property
    def ndim(self):
        return self._ndim_

    @property
    def name(self):
        return self._name_

    def ___PRIVATE_PARSE_type_wrt_metric___(self):
        if self._domain_input_.region_type_wr2_metric is None:
            return TypeWr2MetricGiver('chaotic')(self)
        elif self.name in self._domain_input_.region_type_wr2_metric:
            return TypeWr2MetricGiver(self._domain_input_.region_type_wr2_metric[self.name])(self)
        else:
            return TypeWr2MetricGiver('chaotic')(self)


    @property
    def type_wrt_metric(self):
        return self._type_wrt_metric_

    def ___PRIVATE_set_corner_coordinates___(self, cc):
        """
        Parameters
        ----------
        cc : tuple
            For 2D:
            np.shape(cc) must be (4, 2), and cc[i] corresponds to UL, DL, UR,
            DR corners, and each i is a tuple itself represents (x, y)
            coordinates.

        """
        assert np.shape(cc) == (4, 2), \
            " <Region> : coordinates shape={} wrong.".format(np.shape(cc))
        self._corner_coordinates_ = cc

    @property
    def corner_coordinates(self):
        return self._corner_coordinates_

    @accepts('self', dict)
    def ___PRIVATE_set_edge_types___(self, et):
        """
        Parameters
        ----------
        et : dict
            contain the info of non-straight edges.

        """
        _edge_types_ = [None for _ in range(self.num_edges())]
        for key in et:
            if key.split('-')[0] == self.name:
                _edge_types_[self._edge_name_to_index_(key.split('-')[1])] = et[key]
        self._edge_types_ = tuple(_edge_types_)

    def ___PRIVATE_parse_edge_types___(self):
        """
        Here we get the 4 edge geometries for further getting the region interpolation.

        Attributes
        ----------
        self._edge_geometries_ :
            For 2D: U -> D -> L -> R

        """
        _eg_ = {}
        for i in range(self.num_edges()):
            _et_ = self._edge_types_[i]
            if _et_ is None: # when we do not mention, it is straight, not free.
                _et_ = ('straight',)
            _cc_ = np.array(self.corner_coordinates)[list(self._edge_corner_local_numbering_(i))]
            _eg_[self._edge_index_to_name_(i)] = EdgeGeometryDispatcher(_et_)(_cc_)
        self._edge_geometries_ = _eg_

    @property
    def interpolation(self):
        return self._interpolation_