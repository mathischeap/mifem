# -*- coding: utf-8 -*-



from components.freeze.main import FrozenOnly




class _2dCSCG_Mesh_Elements_DO_FIND(FrozenOnly):
    def __init__(self, elements):
        self._elements_ = elements
        self._freeze_self_()


    def element_at_region_corner(self, region_name, which_corner):
        """

        :param region_name:
        :param which_corner:
        :return:
        """
        mesh = self._elements_._mesh_
        RNS = mesh.domain.regions.names
        assert region_name in RNS, f"region_name = {region_name} is illegal."
        region = mesh.domain.regions[region_name]
        corner_names = region._corner_name_to_index_dict_().keys()
        assert which_corner in corner_names, \
            f"which_corner={which_corner} is illegal!"

        if which_corner[0] in 'LR':
            which_corner = which_corner[::-1]

        id0 = 0 if which_corner[0] == 'U' else -1
        id1 = 0 if which_corner[1] == 'L' else -1

        element_numbering = mesh.___PRIVATE_generate_element_global_numbering___ \
            (number_what=region_name)

        the_corner_element_numbering = element_numbering[id0, id1]

        return the_corner_element_numbering
