# -*- coding: utf-8 -*-

from screws.freeze.main import FrozenOnly
from root.config.main import sIze
import numpy as np


class _2dCSCG_Mesh_DO_FIND(FrozenOnly):
    """A wrapper of all find methods for mesh.do."""
    def __init__(self, meshDO):
        self._DO_ = meshDO
        self._mesh_ = meshDO._mesh_
        self._freeze_self_()

    def region_name_of_element(self, i):
        region_name = None
        for num_elements_accumulation in self._mesh_._num_elements_accumulation_:
            if i < num_elements_accumulation:
                region_name = self._mesh_._num_elements_accumulation_[num_elements_accumulation]
                break
        return region_name

    def slave_of_element(self, i: int) -> int:
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

    def region_name_and_local_indices_of_element(self, i):
        return self._mesh_.___PRIVATE_do_find_region_name_and_local_indices_of_element___(i)

    def reference_origin_and_size_of_element_of_given_local_indices(self, region_name, local_indices):
        origin = [None for _ in range(2)]
        delta = [None for _ in range(2)]
        for i in range(2):
            origin[i] = self._mesh_._element_spacing_[region_name][i][local_indices[i]]
            delta[i] = self._mesh_._element_ratio_[region_name][i][local_indices[i]]
        return tuple(origin), tuple(delta)





    def reference_origin_and_size_of_element(self, i):
        region_name, local_indices = self.region_name_and_local_indices_of_element(i)
        return self.reference_origin_and_size_of_element_of_given_local_indices(
            region_name, local_indices)