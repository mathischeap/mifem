# -*- coding: utf-8 -*-
from components.freeze.base import FrozenOnly



class _2dCSCG_SF_numbering_do_find(FrozenOnly):
    """"""
    def __init__(self, numbering):
        self._numbering_ = numbering
        self._local0SideCache_ = dict()
        self._local_O_1_SideCache_ = dict()
        self._local_I_1_SideCache_ = dict()
        self._freeze_self_()

    def local_dofs_on_element_edge(self, edge_name):
        """

        Parameters
        ----------
        edge_name : str
         {'U', 'D', 'L', 'R'}

        Returns
        -------

        """
        numbering = self._numbering_

        if numbering._sf_.k == 0:
            return self.___PRIVATE_find_0Form_local_dofs_on_element_edge___(edge_name)
        if numbering._sf_.k == 1:
            if numbering._sf_.orientation == 'inner':
                return self.___PRIVATE_find_inner1Form_local_dofs_on_element_edge___(edge_name)
            elif numbering._sf_.orientation == 'outer':
                return self.___PRIVATE_find_outer1Form_local_dofs_on_element_edge___(edge_name)
            else:
                raise Exception
        if numbering._sf_.k == 2:
            raise Exception('volume form has no (local) dofs on element edge.')
        else:
            raise NotImplementedError(f"not coded for {numbering._sf_.k}-form.")

    def ___PRIVATE_find_0Form_local_dofs_on_element_edge___(self, edge_name):
        """"""

        local_numbering = self._numbering_.local

        if edge_name not in self._local0SideCache_:
            if   edge_name == 'U': indices = local_numbering[0][ 0, :]
            elif edge_name == 'D': indices = local_numbering[0][-1, :]
            elif edge_name == 'L': indices = local_numbering[0][:,  0]
            elif edge_name == 'R': indices = local_numbering[0][:, -1]
            else: raise Exception()
            self._local0SideCache_[edge_name] = indices.ravel('F')

        return self._local0SideCache_[edge_name]

    def ___PRIVATE_find_inner1Form_local_dofs_on_element_edge___(self, edge_name):
        """"""
        if edge_name not in self._local_I_1_SideCache_:
            if   edge_name == 'U':
                indices = self._numbering_.local[1][ 0, :]
            elif edge_name == 'D':
                indices = self._numbering_.local[1][-1, :]
            elif edge_name == 'L':
                indices = self._numbering_.local[0][ :, 0]
            elif edge_name == 'R':
                indices = self._numbering_.local[0][ :,-1]
            else:
                raise Exception()
            self._local_I_1_SideCache_[edge_name] = indices.ravel('F')

        return self._local_I_1_SideCache_[edge_name]

    def ___PRIVATE_find_outer1Form_local_dofs_on_element_edge___(self, edge_name):
        """"""
        if edge_name not in self._local_O_1_SideCache_:
            if   edge_name == 'U':
                indices = self._numbering_.local[0][ 0, :]
            elif edge_name == 'D':
                indices = self._numbering_.local[0][-1, :]
            elif edge_name == 'L':
                indices = self._numbering_.local[1][ :, 0]
            elif edge_name == 'R':
                indices = self._numbering_.local[1][ :,-1]
            else:
                raise Exception()
            self._local_O_1_SideCache_[edge_name] = indices.ravel('F')

        return self._local_O_1_SideCache_[edge_name]