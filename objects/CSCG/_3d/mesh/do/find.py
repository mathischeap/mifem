# -*- coding: utf-8 -*-

from root.config.main import sIze, np
from screws.freeze.main import FrozenOnly



class _3dCSCG_Mesh_DO_FIND(FrozenOnly):
    def __init__(self, DO):
        self._DO_ = DO
        self._mesh_ = DO._mesh_
        self._region_names_pool_ = dict()
        self._freeze_self_()

    def region_name_of_element(self, i):
        """Find the regions of ith element.

        Parameters
        ----------
        i

        Returns
        -------

        """
        if i in self._region_names_pool_:
            return self._region_names_pool_[i]
        region_name = None
        for num_elements_accumulation in self._mesh_._num_elements_accumulation_:
            if i < num_elements_accumulation:
                region_name = self._mesh_._num_elements_accumulation_[num_elements_accumulation]
                break
        self._region_names_pool_[i] = region_name
        return region_name

    def region_name_and_local_indices_of_element(self, i):
        """ith global numbered mesh element.

        Parameters
        ----------
        i : int
            The numbering of a global element.

        Returns
        -------

        """
        return self._mesh_.___PRIVATE_do_find_region_name_and_local_indices_of_element___(i)

    def reference_origin_and_size_of_element_of_given_local_indices(self, region_name, local_indices):
        origin = [None for _ in range(self._mesh_.ndim)]
        delta = [None for _ in range(self._mesh_.ndim)]
        for i in range(self._mesh_.ndim):
            origin[i] = self._mesh_._element_spacing_[region_name][i][local_indices[i]]
            delta[i] = self._mesh_._element_ratio_[region_name][i][local_indices[i]]
        return tuple(origin), tuple(delta)

    def reference_origin_and_size_of_element(self, i):
        """
        Find the origin, the UL corner(2D), NWB corner (3D), and the size of
        the ith element in the reference regions [0,1]^ndim.
        """
        region_name, local_indices = self.region_name_and_local_indices_of_element(i)
        return self.reference_origin_and_size_of_element_of_given_local_indices(
            region_name, local_indices)

    def slave_of_element(self, i: int) -> int:
        """Find the core rank of mesh element #i.

        Parameters
        ----------
        i : int
            The number of the mesh element.

        Returns
        -------
        midCore1 : int
            The core the mesh element #`i` is in.

        """
        DISTRI = self._mesh_._element_distribution_
        if isinstance(i, str): i = int(i)
        if sIze <= 6 or not self._mesh_.___is_occupying_all_cores___:
            for nC in range(sIze):
                if i in DISTRI[nC]: return nC
            raise Exception()
        midCore0 = 0
        midCore1 = sIze // 2
        midCore2 = sIze

        while i not in DISTRI[midCore1] and midCore1 - midCore0 > 2 and midCore2 - midCore1 > 2:
            if i > max(DISTRI[midCore1]):
                midCore0 = midCore1
                midCore1 = (midCore0 + midCore2) // 2
            elif i < min(DISTRI[midCore1]):
                midCore2 = midCore1
                midCore1 = (midCore0 + midCore2) // 2
            else:
                raise Exception

        if i in DISTRI[midCore1]:
            return midCore1
        elif i > np.max(DISTRI[midCore1]):
            for noCore in range(midCore1, midCore2):
                if i in DISTRI[noCore]: return noCore
        elif i < np.min(DISTRI[midCore1]):
            for noCore in range(midCore0, midCore1):
                if i in DISTRI[noCore]: return noCore
        else:
            raise Exception

    def element_attach_to_region_side(self, region, side_name):
        """

        :param str region:
        :param str side_name:
        :return:
        """
        EGN1 = self._mesh_.___PRIVATE_generate_element_global_numbering_for_region___(region)
        if side_name == 'N':
            elements = EGN1[ 0, :, :]
        elif side_name == 'S':
            elements = EGN1[-1, :, :]
        elif side_name == 'W':
            elements = EGN1[ :, 0, :]
        elif side_name == 'E':
            elements = EGN1[ :,-1, :]
        elif side_name == 'B':
            elements = EGN1[ :, :, 0]
        elif side_name == 'F':
            elements = EGN1[ :, :,-1]
        else:
            raise Exception()
        return elements