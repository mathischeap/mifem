
from screws.freeze.main import FrozenOnly

from root.config.main import cOmm


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

        pos = cOmm.allgather(pos)

        for _ in pos:
            if _ is not None:
                pos = _
                break

        CHARACTERISTIC_element, side_index = pos
        if CHARACTERISTIC_element in self._elements_._mesh_.elements:
            edge_map = self._elements_._mesh_.edge.elements.map[CHARACTERISTIC_element]
            Region = self._elements_._mesh_.domain.regions.Region
            indices_dict = Region._side_edge_local_numbering_dict_()
            indices = indices_dict[side_index]

            edges = list()
            for _ in indices:
                edges.append(edge_map[_])
        else:
            edges = None

        edges = cOmm.allgather(edges)

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

        MEi = cOmm.allgather(MEi)
        MEj = cOmm.allgather(MEj)

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