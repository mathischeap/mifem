







from screws.frozen import FrozenOnly



class _3dCSCG_Mesh_DO_FIND(FrozenOnly):
    def __init__(self, DO):
        self._DO_ = DO
        self._mesh_ = DO._mesh_
        self._freeze_self_()



    def region_name_of_element(self, i):
        return self._mesh_.___DO_find_region_name_of_element___(i)

    def region_name_and_local_indices_of_element(self, i):
        return self._mesh_.___DO_find_region_name_and_local_indices_of_element___(i)

    def reference_origin_and_size_of_element_of_given_local_indices(self, region_name, local_indices):
        return self._mesh_.___DO_find_reference_origin_and_size_of_element_of_given_local_indices___(
            region_name, local_indices)

    def reference_origin_and_size_of_element(self, i):
        return self._mesh_.___DO_find_reference_origin_and_size_of_element___(i)

    def slave_of_element(self, i):
        """Find the core rank of mesh element #i."""
        return self._mesh_.___DO_find_slave_of_element___(i)

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
