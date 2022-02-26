

from screws.frozen import FrozenOnly
from _3dCSCG.mesh.DO.FIND import _3dCSCG_Mesh_DO_FIND


class _3dCSCG_Mesh_DO(FrozenOnly):
    def __init__(self, mesh):
        self._mesh_ = mesh
        self._FIND_ = _3dCSCG_Mesh_DO_FIND(self)
        self._freeze_self_()

    def RESET_cache(self):
        self._mesh_.RESET_cache()

    def parse_element_side_pair(self, eP):
        return self._mesh_.___DO_parse_element_side_pair___(eP)



    def FIND_region_name_of_element(self, i):
        return self._mesh_.___DO_find_region_name_of_element___(i)

    def FIND_region_name_and_local_indices_of_element(self, i):
        return self._mesh_.___DO_find_region_name_and_local_indices_of_element___(i)

    def FIND_reference_origin_and_size_of_element_of_given_local_indices(self, region_name, local_indices):
        return self._mesh_.___DO_find_reference_origin_and_size_of_element_of_given_local_indices___(
            region_name, local_indices)

    def FIND_reference_origin_and_size_of_element(self, i):
        return self._mesh_.___DO_find_reference_origin_and_size_of_element___(i)

    def FIND_slave_of_element(self, i):
        """Find the core rank of mesh element #i."""
        return self._mesh_.___DO_find_slave_of_element___(i)

    def FIND_element_attach_to_region_side(self, region, side_name):
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

    @property
    def FIND(self):
        return self._FIND_


    def regionwsie_stack(self, *args):
        return self._mesh_.___DO_regionwsie_stack___(*args)

