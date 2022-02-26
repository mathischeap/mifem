# -*- coding: utf-8 -*-


import numpy as np
from root.config import cOmm
from screws.frozen import FrozenOnly
from _3dCSCG.mesh.regions.main import Regions
from _3dCSCG.mesh.regions.region.main import Region
from _3dCSCG.mesh.domain.visualize.main import _3dCSCG_Domain_Visualize
from _3dCSCG.mesh.domain.boundaries import _3dCSCG_Boundaries


class _3dCSCG_Domain(FrozenOnly):
    """We have the whole ``_3dCSCG_Domain`` (all same) in all cores. This
    is very important.
    """
    def __init__(self, di):
        self.___PRIVATE_parse_domain_input___(di)
        self.___PRIVATE_parse_interpolators___()
        self.___PRIVATE_generate_regions___()
        self.___PRIVATE_generate_global_corner_numbering___()
        self.___PRIVATE_generate_region_map___()
        v = 0
        for rn in self.regions.names:
            v += self.regions(rn).volume
        self._volume_ = v
        self._visualize_ = None
        self._boundaries_ = None
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

    @property
    def domain_input(self):
        return self._domain_input_

    @property
    def ndim(self):
        return self.domain_input.ndim

    @property
    def volume(self):
        return self._volume_

    @property
    def visualize(self):
        if self._visualize_ is None:
            self._visualize_ = _3dCSCG_Domain_Visualize(self)
        return self._visualize_

    @property
    def boundaries(self):
        if self._boundaries_ is None:
            self._boundaries_ = _3dCSCG_Boundaries(self)
        return self._boundaries_

    @property
    def interpolators(self):
        """
        Returns
        -------
        self._interpolator_ : dict

        """
        return self._interpolators_

    @property
    def name(self):
        return self.domain_input.domain_name

    def __len__(self):
        """
        The length of a domain, len(domain), is defined to be the number of
        regions the domain has.

        """
        return self._num_regions_

    @property
    def regions(self):
        """
        Returns
        -------
        self._regions_ : dict
            A dict that contains all regions. Keys: regions' names. Values: the
            regions instances.

        """
        return self._regions_



    def ___PRIVATE_parse_domain_input___(self, di):
        assert di.ndim == 3, " <Domain> : I need 3d DomainInput."
        self._domain_input_ = di
        self._num_regions_ = len(di.region_corner_coordinates)
        self._num_boundaries_ = len(di.boundary_region_sides)

        # important: we have an option to customize the regions sequence which is important for like element numbering.
        if di.region_sequence is None:
            self._region_names_ = tuple(di.region_corner_coordinates.keys())
        else:
            self._region_names_ = di.region_sequence

        self._boundary_names_ = tuple(di.boundary_region_sides.keys())

    def ___PRIVATE_parse_interpolators___(self):
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
        if isinstance(self.domain_input.region_interpolators, str):
            self._interpolators_ = dict()
            for rn in self._region_names_:
                self._interpolators_[rn] = self.domain_input.region_interpolators
        elif isinstance(self.domain_input.region_interpolators, dict):
            self._interpolators_ = self.domain_input.region_interpolators
        else:
            raise Exception(" <Domain> ")
        assert set(self._interpolators_.keys()) == set(self._region_names_), \
            " <Domain> : I need interpolator for every regions. "

    def ___PRIVATE_generate_regions___(self):
        """
        We use this property to parse self.domain_input. Then we can get
        several regions (instances of Region2D or Region3D).

        """
        cOmm.barrier()
        _regions_ = dict()
        for rn in self._region_names_:
            _R_ = Region(
                self.ndim,
                rn,
                self.domain_input.region_corner_coordinates[rn],
                self.domain_input.region_side_types,
                self.interpolators[rn],
                self.domain_input)
            _regions_[rn] = _R_
        self._regions_ = Regions(self, _regions_)

    def ___PRIVATE_generate_global_corner_numbering___(self):
        """
        Here we parse our regions and try to get the global_corner_numbering
        and something else.

        `_region_corner_coordinates_pool_` and `_global_region_corner_numbering_` are
        not very useful because when two sides at the same location but are on domain
        boundary (like a crack), the global corner numbering will be wrong. So,
        generally, we do not use it.

        Attributes
        ----------
        self._region_corner_coordinates_pool_ : tuple
            A tuple where self._region_corner_pool_[i] indicates the coordinates of the
            regions corner globally numbered as i.
        self._global_region_corner_numbering_ : dict
            A dict whose names indicate regions, and whose valune[j].ravel('F')
            indicates the golobal numbering of [NWB, SWB, NEB, SEB, NWF, SWF,
            NEB, SEB] [j] corner.

        """
        _rcgn_ = dict()
        _corner_pool_ = list()
        _current_number_ = 0
        for rn in self.regions():  # go through all regions
            _rcgn_[rn] = ()
            for i in range(self.regions(rn).num_corners()):  # go through all corners.
                _corner_ = self.regions(rn).corner_coordinates[i]
                _in_, _index_ = self.___PRIVATE_check_duplicated_corner___(_corner_pool_, _corner_)
                if _in_:
                    _rcgn_[rn] += (_index_,)
                else:
                    _rcgn_[rn] += (_current_number_,)
                    _current_number_ += 1
                    _corner_pool_.append(_corner_)
            _rcgn_[rn] = np.array(_rcgn_[rn]).reshape((2, 2, 2), order='F')
        self._region_corner_coordinates_pool_ = tuple(_corner_pool_)
        self._region_corner_global_numbering_ = _rcgn_

    def ___PRIVATE_generate_region_map___(self):
        """
        Generate the regions map and something else.

        Attributes
        ----------
        self._region_internal_side_pairs_ : tuple
            A tuple of sets of two elements. Each set represents two internal
            regions sides which are paired up.
        self._region_sides_on_domain_boundaries_ : dict
        self._region_map_ : dict

        """
        _rm_ = dict()
        # we first find the internal pairing___________________________________
        rn = None

        for rn in self.regions():  # go through all regions
            _rm_[rn] = [list() for _ in range(self.regions(rn).num_sides())]
            for i in range(self.regions(rn).num_sides()):  # go through all 6 sides of each regions.
                self_corner_indices = self.___PRIVATE_found_side_corner_global_numbering___(rn, i)
                for rnrn in self.regions().keys():  # go through all regions except self
                    if rnrn != rn:
                        for ii in range(self.regions(rn).num_sides()):  # go through all 6 sides of the regions.
                            other_corner_indices = self.___PRIVATE_found_side_corner_global_numbering___(rnrn, ii)
                            if np.all(self_corner_indices == other_corner_indices):
                                _rm_[rn][i].append(rnrn)

        # We then find the sides on the domain boundaries______________________
        for bn in self._boundary_names_:  # we of course go through all boundaries
            for db_i in self.domain_input.boundary_region_sides[bn]:
                _region_name_, _region_side_ = db_i.split('-')
                _rm_[_region_name_][self.regions(rn)._side_name_to_index_(_region_side_)].append(bn)

        # Now we check the regions map and extract more info____________________
        # We first check each side only appears at one place and extract the sides on domain boundaries
        _rsodb_ = {}
        for rn in self._region_names_:  # go through all regions
            _rsodb_[rn] = [0 for _ in range(self.regions(rn).num_sides())]
            for i in range(self.regions(rn).num_sides()):  # go through all 6 sides of each regions.
                try:
                    assert np.shape(_rm_[rn][i]) == (1,)
                except AssertionError:
                    assert np.shape(_rm_[rn][i]) == (2,), \
                        " <Domain3D> : region_map[{}][{}] = {} is wrong, check the domain_input.".format(
                            rn, i, _rm_[rn][i])
                    rmrni0, rmrni1 = _rm_[rn][i]
                    if rmrni0 in self._region_names_:
                        assert rmrni1 in self._boundary_names_
                        rmrnib = rmrni1
                        rmrnir = rmrni0
                    elif rmrni0 in self._boundary_names_:
                        assert rmrni1 in self._region_names_
                        rmrnib = rmrni0
                        rmrnir = rmrni1
                    else:
                        raise Exception()
                    assert np.shape(_rm_[rmrnir][{0: 1, 1: 0, 2: 3, 3: 2, 4: 5, 5: 4}[i]]) == (2,)
                    assert rn in _rm_[rmrnir][{0: 1, 1: 0, 2: 3, 3: 2, 4: 5, 5: 4}[i]]
                    assert rmrnib in _rm_[rmrnir][{0: 1, 1: 0, 2: 3, 3: 2, 4: 5, 5: 4}[i]]
                    _rm_[rn][i] = [rmrnib, ]
                    _rm_[rmrnir][{0: 1, 1: 0, 2: 3, 3: 2, 4: 5, 5: 4}[i]] = [rmrnib, ]
                _rm_[rn][i] = _rm_[rn][i][0]
                if _rm_[rn][i] in self._boundary_names_:
                    _rsodb_[rn][i] = 1
                    assert _rm_[rn][i] not in self._region_names_, \
                        " <Domain3D> : region_map[{}][{}] = {} is wrong, check the domain_input.".format(
                            rn, i, _rm_[rn][i])
                else:
                    assert _rm_[rn][i] in self._region_names_, \
                        " <Domain3D> : region_map[{}][{}] = {} is wrong, check the domain_input.".format(
                            rn, i, _rm_[rn][i])

            _rsodb_[rn] = tuple(_rsodb_[rn])
            _rm_[rn] = tuple(_rm_[rn])

        # now we check internal sides are correctly paired up and extract the pairing info
        _region_side_correct_pairing_ = ({0, 1}, {2, 3}, {4, 5})
        _risp_ = ()
        for rn in self._region_names_:  # go through all regions
            for i in range(self.regions(rn).num_sides()):  # go through all sides of each regions.
                if not _rsodb_[rn][i]:
                    assert {i, _rm_[_rm_[rn][i]].index(rn)} in _region_side_correct_pairing_, \
                        " <Domain3D> : regions['{}']-side[{}] is paired up with regions['{}']-side[{}]".format(
                            rn, i, _rm_[rn][i], _rm_[_rm_[rn][i]].index(rn))
                    current_pair = {
                        rn + '-' + self.regions(rn)._side_index_to_name_(i),
                        str(_rm_[rn][i]) + '-' + self.regions(rn)._side_index_to_name_(_rm_[_rm_[rn][i]].index(rn))}
                    if current_pair not in _risp_:
                        _risp_ += (current_pair,)

        self._region_internal_side_pairs_ = _risp_
        self._region_sides_on_domain_boundaries_ = _rsodb_
        self._region_map_ = _rm_

    def ___PRIVATE_found_side_corner_global_numbering___(self, region_name, side_index):
        """
        Given a tuple(size 2 or 4) of 2 or 4 local numberings of 2 or 4 corners
        of a side (corner_local_numbering) of a regions (region_name), we return
        the 2 or 4  global numberings.

        Parameters
        ----------
        region_name : str
        side_index : int

        Returns
        -------
        output : ndarray

        """
        return self._region_corner_global_numbering_[region_name].ravel('F')[
            list(self.regions(region_name)._side_corner_local_numbering_(side_index))]

    @staticmethod
    def ___PRIVATE_check_duplicated_corner___(pool, corner):
        """
        Here we check if a corner: corner is already recorded in pool.

        Parameters
        ----------
        pool : list
        corner : tuple

        Returns
        -------
        _in_ : pool
        _index_ : int

        """
        if corner in pool:
            return True, pool.index(corner)
        for i in pool:
            if np.sqrt((i[0] - corner[0]) ** 2 + (i[1] - corner[1]) ** 2 + (i[2] - corner[2]) ** 2) <= 1e-13:
                return True, i
        return False, -1



    @property
    def IS_periodic(self):
        return True if len(self.domain_input.periodic_boundary_pairs) != 0 else False