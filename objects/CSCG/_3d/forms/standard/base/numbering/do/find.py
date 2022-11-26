# -*- coding: utf-8 -*-
from components.freeze.main import FrozenOnly


class _3dCSCG_Standard_Form_Numbering_DO_FIND_(FrozenOnly):
    """"""

    def __init__(self, numbering_do):
        """"""
        self._numbering_ = numbering_do._numbering_
        self._local0SideCache_ = dict()
        self._local1SideCache_ = dict()
        self._local2SideCache_ = dict()
        self._local_1_ECE_Cache_ = dict()
        self._local_0_ECE_Cache_ = dict()
        self._freeze_self_()

    def dofs_on_element_side(self, element, side_name, GM=None):
        """

        :param element:
        :param side_name:
        :param GM: If GM is None, we use self.GM, otherwise, we use this given GM.
        :return:
        """
        numbering = self._numbering_

        if GM is None:
            GM = numbering.gathering
        if numbering._sf_.k == 0:
            return numbering.___PRIVATE_find_0Form_dofs_on_element_side___(element, side_name, GM)
        if numbering._sf_.k == 1:
            return numbering.___PRIVATE_find_1Form_dofs_on_element_side___(element, side_name, GM)
        if numbering._sf_.k == 2:
            return numbering.___PRIVATE_find_2Form_dofs_on_element_side___(element, side_name, GM)
        if numbering._sf_.k == 3:
            raise Exception('volume form has no dofs on element side.')
        else:
            raise NotImplementedError(f"not coded for {numbering._sf_.k}-form.")

    def local_dofs_on_element_side(self, side_name):
        """

        :param side_name:
        :return:
        """
        numbering = self._numbering_

        if numbering._sf_.k == 0:
            return self.___PRIVATE_find_0Form_local_dofs_on_element_side___(side_name)
        if numbering._sf_.k == 1:
            return self.___PRIVATE_find_1Form_local_dofs_on_element_side___(side_name)
        if numbering._sf_.k == 2:
            return self.___PRIVATE_find_2Form_local_dofs_on_element_side___(side_name)
        if numbering._sf_.k == 3:
            raise Exception('volume form has no (local) dofs on element side.')
        else:
            raise NotImplementedError(f"not coded for {numbering._sf_.k}-form.")

    def ___PRIVATE_find_0Form_local_dofs_on_element_side___(self, side_name):
        """"""
        local_numbering = self._numbering_.local[0]

        if side_name not in self._local0SideCache_:
            if   side_name == 'N': indices = local_numbering[ 0, :, :]
            elif side_name == 'S': indices = local_numbering[-1, :, :]
            elif side_name == 'W': indices = local_numbering[ :, 0, :]
            elif side_name == 'E': indices = local_numbering[ :,-1, :]
            elif side_name == 'B': indices = local_numbering[ :, :, 0]
            elif side_name == 'F': indices = local_numbering[ :, :,-1]
            else: raise Exception()
            self._local0SideCache_[side_name] = indices.ravel('F')

        return self._local0SideCache_[side_name]

    def ___PRIVATE_find_1Form_local_dofs_on_element_side___(self, side_name):
        """"""
        if side_name not in self._local1SideCache_:
            if   side_name == 'N':
                indices1 = self._numbering_.local[1][ 0, :, :]
                indices2 = self._numbering_.local[2][ 0, :, :]
            elif side_name == 'S':
                indices1 = self._numbering_.local[1][-1, :, :]
                indices2 = self._numbering_.local[2][-1, :, :]
            elif side_name == 'W':
                indices1 = self._numbering_.local[0][ :, 0, :]
                indices2 = self._numbering_.local[2][ :, 0, :]
            elif side_name == 'E':
                indices1 = self._numbering_.local[0][ :,-1, :]
                indices2 = self._numbering_.local[2][ :,-1, :]
            elif side_name == 'B':
                indices1 = self._numbering_.local[0][ :, :, 0]
                indices2 = self._numbering_.local[1][ :, :, 0]
            elif side_name == 'F':
                indices1 = self._numbering_.local[0][ :, :,-1]
                indices2 = self._numbering_.local[1][ :, :,-1]
            else: raise Exception()
            self._local1SideCache_[side_name] = list()
            self._local1SideCache_[side_name].extend(indices1.ravel('F'))
            self._local1SideCache_[side_name].extend(indices2.ravel('F'))

        return self._local1SideCache_[side_name]

    def ___PRIVATE_find_2Form_local_dofs_on_element_side___(self, side_name):
        """"""
        if side_name not in self._local2SideCache_:

            if   side_name == 'N': indices = self._numbering_.local[0][ 0, :, :]
            elif side_name == 'S': indices = self._numbering_.local[0][-1, :, :]
            elif side_name == 'W': indices = self._numbering_.local[1][ :, 0, :]
            elif side_name == 'E': indices = self._numbering_.local[1][ :,-1, :]
            elif side_name == 'B': indices = self._numbering_.local[2][ :, :, 0]
            elif side_name == 'F': indices = self._numbering_.local[2][ :, :,-1]
            else: raise Exception()
            self._local2SideCache_[side_name] = indices.ravel('F')

        return self._local2SideCache_[side_name]

    def local_dofs_on_element_corner_edge(self, corner_edge_name):
        """"""
        names = ['WB', 'EB', 'WF', 'EF', 'NB', 'SB', 'NF', 'SF', 'NW', 'SW', 'NE', 'SE']
        assert corner_edge_name in names, f"corner_edge_name={corner_edge_name} is illegal."
        k = self._numbering_._sf_.k
        if k == 0:
            return self.___PRIVATE_find_S0F_local_dofs_on_element_corner_edge___(corner_edge_name)
        elif k == 1:
            return self.___PRIVATE_find_S1F_local_dofs_on_element_corner_edge___(corner_edge_name)
        else:
            raise Exception()

    def ___PRIVATE_find_S0F_local_dofs_on_element_corner_edge___(self, corner_edge_name):
        """"""
        if corner_edge_name not in self._local_0_ECE_Cache_:

            if corner_edge_name == 'WB':
                indices = self._numbering_.local[0][:, 0, 0]
            elif corner_edge_name == 'EB':
                indices = self._numbering_.local[0][:, -1, 0]
            elif corner_edge_name == 'WF':
                indices = self._numbering_.local[0][:, 0, -1]
            elif corner_edge_name == 'EF':
                indices = self._numbering_.local[0][:, -1, -1]
            elif corner_edge_name == 'NB':
                indices = self._numbering_.local[0][0, :, 0]
            elif corner_edge_name == 'SB':
                indices = self._numbering_.local[0][-1, :, 0]
            elif corner_edge_name == 'NF':
                indices = self._numbering_.local[0][0, :, -1]
            elif corner_edge_name == 'SF':
                indices = self._numbering_.local[0][-1, :, -1]
            elif corner_edge_name == 'NW':
                indices = self._numbering_.local[0][0, 0, :]
            elif corner_edge_name == 'SW':
                indices = self._numbering_.local[0][-1, 0, :]
            elif corner_edge_name == 'NE':
                indices = self._numbering_.local[0][0, -1, :]
            elif corner_edge_name == 'SE':
                indices = self._numbering_.local[0][-1, -1, :]
            else:
                raise Exception()

            self._local_0_ECE_Cache_[corner_edge_name] = indices

        return self._local_0_ECE_Cache_[corner_edge_name]

    def ___PRIVATE_find_S1F_local_dofs_on_element_corner_edge___(self, corner_edge_name):
        """"""
        if corner_edge_name not in self._local_1_ECE_Cache_:

            if corner_edge_name == 'WB':
                indices = self._numbering_.local[0][:, 0, 0]
            elif corner_edge_name == 'EB':
                indices = self._numbering_.local[0][:, -1, 0]
            elif corner_edge_name == 'WF':
                indices = self._numbering_.local[0][:, 0, -1]
            elif corner_edge_name == 'EF':
                indices = self._numbering_.local[0][:, -1, -1]
            elif corner_edge_name == 'NB':
                indices = self._numbering_.local[1][0, :, 0]
            elif corner_edge_name == 'SB':
                indices = self._numbering_.local[1][-1, :, 0]
            elif corner_edge_name == 'NF':
                indices = self._numbering_.local[1][0, :, -1]
            elif corner_edge_name == 'SF':
                indices = self._numbering_.local[1][-1, :, -1]
            elif corner_edge_name == 'NW':
                indices = self._numbering_.local[2][0, 0, :]
            elif corner_edge_name == 'SW':
                indices = self._numbering_.local[2][-1, 0, :]
            elif corner_edge_name == 'NE':
                indices = self._numbering_.local[2][0, -1, :]
            elif corner_edge_name == 'SE':
                indices = self._numbering_.local[2][-1, -1, :]
            else:
                raise Exception()

            self._local_1_ECE_Cache_[corner_edge_name] = indices

        return self._local_1_ECE_Cache_[corner_edge_name]