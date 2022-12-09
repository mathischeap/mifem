# -*- coding: utf-8 -*-
from components.freeze.main import FrozenOnly

from root.config.main import COMM

class _3dCSCG_Trace_Elements_Do_Find(FrozenOnly):
    """We find some specific groups of elements."""

    def __init__(self, trace_elements):
        self._elements_ = trace_elements
        self._freeze_self_()

    def element_between_two_mesh_elements(self, i, j):
        """We try to find the trace element between mesh elements #`i` and #`j`.

        Parameters
        ----------
        i : int
        j : int

        Returns
        -------

        """
        return self._elements_._mesh_.elements.do.find.trace_element_between_two_elements(i, j)

    def edge_elements_surrounding_element(self, i):
        """We try to find the four edge-elements surrounding trace element #`i`.

        Parameters
        ----------
        i : int
            The trace-element.

        Returns
        -------

        """

        if i in self._elements_:
            CHARACTERISTIC_element = self._elements_[i].CHARACTERISTIC_element
            MAP = self._elements_.map[CHARACTERISTIC_element]
            side_index = MAP.index(i)
            pos = (CHARACTERISTIC_element, side_index)
        else:
            pos = None

        pos = COMM.allgather(pos)

        for _ in pos:
            if _ is not None:
                pos = _
                break

        edge_elements = self._elements_._mesh_.edge.elements # in case the edge elements were not made.
        CHARACTERISTIC_element, side_index = pos
        if CHARACTERISTIC_element in self._elements_._mesh_.elements:
            edge_map = edge_elements.map[CHARACTERISTIC_element]
            Region = self._elements_._mesh_.domain.regions.Region
            indices_dict = Region._side_edge_local_numbering_dict_()
            indices = indices_dict[side_index]

            edges = list()
            for _ in indices:
                edges.append(edge_map[_])
        else:
            edges = None

        edges = COMM.allgather(edges)

        for _ in edges:
            if _ is not None:
                edges = _
                break

        return edges

    def elements_attached_to_edge_element(self, i):
        """Find the at most 4 trace elements that are attached to the edge element #`i`.

        Parameters
        ----------
        i : int

        Returns
        -------

        """
        return self._elements_._mesh_.edge.elements.do.find.trace_elements_attached_to_element(i)

    def elements_attached_to_node_element(self, i):
        """Find the trace elements that are attached to a node element #`i`.

        For example, for a typical internal node element, there will be 12 trace elements attached
        to it.

        Parameters
        ----------
        i : int
            The node element #`i`.

        Returns
        -------

        """
        return self._elements_._mesh_.node.elements.do.find.trace_elements_attached_to_element(i)

    def node_elements_of_element(self, i):
        """Return a list of the (four) node elements of trace element #i.

        Parameters
        ----------
        i

        Returns
        -------

        """

        if i in self._elements_:
            DICT = self._elements_._mesh_.domain.regions.Region._side_corner_local_numbering_dict_()

            te = self._elements_[i]
            CHARACTERISTIC_element = te.CHARACTERISTIC_element
            CHARACTERISTIC_side = te.CHARACTERISTIC_side

            NODE_MAP = self._elements_._mesh_.node.elements.map[CHARACTERISTIC_element]

            corner_indices = DICT['NSWEBF'.index(CHARACTERISTIC_side)]

            node_elements = list()
            for ci in corner_indices:
                node_elements.append(NODE_MAP[ci])
        else:
            node_elements = None

        node_elements = COMM.allgather(node_elements)

        NE = None
        for _ in node_elements:
            if _ is not None:
                if NE is None:
                    NE = _
                else:
                    assert NE == _

        return NE

    def mesh_element_shared_by_elements(self, i, j):
        """Find the mesh element shared by the two trace elements #i, #j. And if they do not share
        a mesh element, return None.

        Parameters
        ----------
        i
        j

        Returns
        -------

        """
        if i in self._elements_:
            TEi = self._elements_[i]
            if TEi.NON_CHARACTERISTIC_position[0].isnumeric():
                MEi = [int(TEi.CHARACTERISTIC_position[:-1]), int(TEi.NON_CHARACTERISTIC_position[:-1])]
            else:
                MEi = [int(TEi.CHARACTERISTIC_position[:-1]),]
        else:
            MEi = None

        if j in self._elements_:
            TEj = self._elements_[j]
            if TEj.NON_CHARACTERISTIC_position[0].isnumeric():
                MEj = [int(TEj.CHARACTERISTIC_position[:-1]), int(TEj.NON_CHARACTERISTIC_position[:-1])]
            else:
                MEj = [int(TEj.CHARACTERISTIC_position[:-1]),]
        else:
            MEj = None

        MEi = COMM.allgather(MEi)
        MEj = COMM.allgather(MEj)

        for _ in MEi:
            if _ is not None:
                MEi = _
                break

        for _ in MEj:
            if _ is not None:
                MEj = _
                break

        for ei in MEi:
            if ei in MEj:
                return ei # we find one shared mesh element, return!

        return None # return None as we do not find a shared mesh element