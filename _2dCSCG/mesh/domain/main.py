# -*- coding: utf-8 -*-
"""
INTRO

@author: Yi Zhang. Created on Tue May 21 14:10:36 2019
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft
         Delft, Netherlands

"""
from SCREWS.frozen import FrozenOnly
from _2dCSCG.mesh.region.main import Region, Regions
from _2dCSCG.mesh.domain.visualize import _2dCSCG_Domain_Visualize
from _2dCSCG.mesh.domain.boundaries import _2dCSCG_Domain_Boundaries


class _2dCSCG_Domain(FrozenOnly):
    def __init__(self, di):
        """
        Parameters
        ---------
        di : DomainInput
            The DomainInput instance.

        """
        self.___PRIVATE_parse_domain_input___(di)
        self.___PRIVATE_parse_interpolators_()
        self.___PRIVATE_generate_regions___()
        self._visualize_ = _2dCSCG_Domain_Visualize(self) # will only do thing in master core.
        self._boundaries_ = _2dCSCG_Domain_Boundaries(self)
        self.___define_parameters___ = None
        self._freeze_self_()

    @property
    def parameters(self):
        return self.___define_parameters___

    def __eq__(self, other):
        for key in self.parameters:
            if key in self.domain_input.internal_parameters:
                pass
            else:
                if self.parameters[key] != other.parameters[key]:
                    return False
        return True


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
        Here only get the interpolator names actually. The interpolation class
        will be obtained in the region itself.

        We get the `_interpolators_` from the `domain_input`:
        `self.domain_input.region_interpolators`. If
        `self.domain_input.region_interpolators` is str, then all regions use
        will use the same interpolator named after this str. If
        `self.domain_input.region_interpolators` is dict, then we take this
        dict directly. Else, we raise Exception.

        Attributes
        ----------
        self._interpolators_ : dict
            The dict whose keys are region names and values are the
            interpolator names the interpolators of regions.

        """
        if isinstance(self.domain_input.region_interpolators, str):
            self._interpolators_ = {}
            for rn in self._region_names_:
                self._interpolators_[rn] = self.domain_input.region_interpolators
        elif isinstance(self.domain_input.region_interpolators, dict):
            self._interpolators_ = self.domain_input.region_interpolators
        else:
            raise Exception(" <Domain> ")
        assert set(self._interpolators_.keys()) == set(self._region_names_), \
            " <Domain> : I need interpolator for every region. "

    def ___PRIVATE_generate_regions___(self):
        """
        We use this property to parse self.domain_input. Then we can get
        several regions (instances of Region2D or Region3D).

        """
        _regions_ = {}
        for rn in self._region_names_:
            _regions_[rn] = self.Region(
                self.ndim,
                rn,
                self.domain_input.region_corner_coordinates[rn],
                self.domain_input.region_edge_types,
                self.interpolators[rn],
                self.domain_input)
        self._regions_ = Regions(self, _regions_)

    @property
    def ndim(self):
        return self.domain_input.ndim

    @property
    def domain_input(self):
        return self._domain_input_

    @property
    def Region(self):
        return Region

    @property
    def name(self):
        return self.domain_input.domain_name

    @property
    def visualize(self):
        return self._visualize_

    @property
    def boundaries(self):
        return self._boundaries_


    def __len__(self):
        """
        The length of a domain, len(domain), is defined to be the number of
        regions the domain has.

        """
        return self._num_regions_


    @property
    def interpolators(self):
        """
        Returns
        -------
        self._interpolator_ : dict

        """
        return self._interpolators_

    @property
    def regions(self):
        """
        Returns
        -------
        self._regions_ : dict
            A dict that contains all regions. Keys: regions' names. Values: the
            region instances.

        """
        return self._regions_


    @property
    def IS_periodic(self):
        return True if len(self.domain_input.periodic_boundary_pairs) != 0 else False