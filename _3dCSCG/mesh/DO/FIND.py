







from screws.frozen import FrozenOnly



class _3dCSCG_Mesh_DO_FIND(FrozenOnly):
    def __init__(self, DO):
        self._DO_ = DO
        self._freeze_self_()



    def region_name_of_element(self, i):
        return self._DO_.FIND_region_name_of_element(i)

    def region_name_and_local_indices_of_element(self, i):
        return self._DO_.FIND_region_name_and_local_indices_of_element(i)

    def reference_origin_and_size_of_element_of_given_local_indices(self, region_name, local_indices):
        return self._DO_.FIND_reference_origin_and_size_of_element_of_given_local_indices(
            region_name, local_indices)

    def reference_origin_and_size_of_element(self, i):
        return self._DO_.FIND_reference_origin_and_size_of_element(i)

    def slave_of_element(self, i):
        """Find the core rank of mesh element #i."""
        return self._DO_.FIND_slave_of_element(i)

    def element_attach_to_region_side(self, region, side_name):
        """

        :param str region:
        :param str side_name:
        :return:
        """
        return self._DO_.FIND_element_attach_to_region_side(region, side_name)
