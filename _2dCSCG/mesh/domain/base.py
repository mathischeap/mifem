
from screws.freeze.main import FrozenOnly
from _2dCSCG.mesh.domain.regions.main import Regions
from _2dCSCG.mesh.domain.regions.region.main import Region


class _2dCSCG_DomainBase(FrozenOnly):
    def __init__(self, di):
        self.___PRIVATE_parse_domain_input___(di)
        self.___PRIVATE_parse_interpolators_()
        self.___PRIVATE_generate_regions___()


    def ___PRIVATE_parse_domain_input___(self, di):
        assert di.ndim == 2, " <Domain> : I need 2d DomainInput."
        self._domain_input_ = di
        self._num_boundaries_ = len(di.boundary_region_edges)
        self._boundary_names_ = tuple(di.boundary_region_edges.keys())
        self._num_regions_ = len(di.region_corner_coordinates)
        if di.region_sequence is None:
            self._region_names_ = tuple(di.region_corner_coordinates.keys())
        else:
            self._region_names_ = di.region_sequence

    def ___PRIVATE_parse_interpolators_(self):
        """
        Here only get the interpolator names. The interpolation class
        will be obtained in the regions itself.

        We get the `_interpolators_` from the `domain_input`:
        `self.domain_input.region_interpolators`. If
        `self.domain_input.region_interpolators` is str, then all regions use
        will use the same interpolator named after this str. If
        `self.domain_input.region_interpolators` is dict, then we take this
        dict directly. Else, we raise Exception.

        Attributes
        ----------
        self._interpolators_ : dict
            The dict whose keys are regions names and values are the
            interpolator names the interpolators of regions.

        """
        if isinstance(self._domain_input_.region_interpolators, str):
            self._interpolators_ = {}
            for rn in self._region_names_:
                self._interpolators_[rn] = self._domain_input_.region_interpolators
        elif isinstance(self._domain_input_.region_interpolators, dict):
            self._interpolators_ = self._domain_input_.region_interpolators
        else:
            raise Exception(" <Domain> ")
        assert set(self._interpolators_.keys()) == set(self._region_names_), \
            " <Domain> : I need interpolator for every regions. "

    def ___PRIVATE_generate_regions___(self):
        """
        We use this property to parse self.domain_input. Then we can get
        several regions (instances of Region2D or Region3D).

        """
        _regions_ = {}
        for rn in self._region_names_:
            _regions_[rn] = Region(
                2, #ndim
                rn,
                self._domain_input_.region_corner_coordinates[rn],
                self._domain_input_.region_edge_types,
                self._interpolators_[rn],
                self._domain_input_)
        self._regions_ = Regions(self, _regions_)