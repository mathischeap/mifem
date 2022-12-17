# -*- coding: utf-8 -*-
""""""
from components.freeze.main import FrozenOnly
from root.config.main import COMM


class _3dCSCG_Mesh_Elements_DO_FIND(FrozenOnly):
    def __init__(self, elements):
        self._elements_ = elements
        self._freeze_self_()

    def slave_of_element(self, i):
        """Find the core rank of mesh element #i."""
        return self._elements_._mesh_.do.find.slave_of_element(i)

    def element_at_region_corner(self, region_name, which_corner):
        """Return the mesh element numbering at the given corner of the domain region.

        Parameters
        ----------
        region_name : str
            The region name
        which_corner : str
            The corner name, for example, `'NWB'`, `'SEF'`.

        Returns
        -------
        the_corner_element_numbering : int
            The numbering of the mesh element at the region corner (No matter if this mesh
            element is in this core).

        """
        mesh = self._elements_._mesh_
        RNS = mesh.domain.regions.names
        assert region_name in RNS, f"region_name = {region_name} is illegal."
        region = mesh.domain.regions[region_name]
        corner_names = region._corner_name_to_index_dict_().keys()
        assert which_corner in corner_names, \
            f"which_corner={which_corner} is illegal!"

        if 'N' in corner_names:
            id0 = 0
        else:
            assert 'S' in corner_names
            id0 = -1

        if 'W' in corner_names:
            id1 = 0
        else:
            assert 'E' in corner_names
            id1 = -1

        if 'B' in corner_names:
            id2 = 0
        else:
            assert 'F' in corner_names
            id2 = -1

        element_numbering = mesh.___PRIVATE_generate_element_global_numbering___(number_what=region_name)

        the_corner_element_numbering = element_numbering[id0, id1, id2]

        return the_corner_element_numbering

    def trace_element_between_two_elements(self, i, j):
        """We try to find the trace element between mesh elements #`i` and #`j`.

        Parameters
        ----------
        i : int
        j : int

        Returns
        -------

        """
        if i in self._elements_:
            TEsi = self._elements_._mesh_.trace.elements.map[i]
        else:
            TEsi = None

        if j in self._elements_:
            TEsj = self._elements_._mesh_.trace.elements.map[j]
        else:
            TEsj = None

        TEsi = COMM.allgather(TEsi)
        TEsj = COMM.allgather(TEsj)

        for _ in TEsi:
            if _ is not None:
                TEsi = _
                break

        for _ in TEsj:
            if _ is not None:
                TEsj = _
                break

        T = list()
        for t in TEsi:
            if t in TEsj:
                T.append(t)

        return T

    def element_shared_by_trace_elements(self, i, j):
        """

        Parameters
        ----------
        i
        j

        Returns
        -------

        """
        return self._elements_._mesh_.trace.elements.do.find.mesh_element_shared_by_elements(i, j)
