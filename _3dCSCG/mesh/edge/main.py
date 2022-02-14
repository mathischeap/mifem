# -*- coding: utf-8 -*-
"""
The edge elements of a mesh.
"""


import sys
if './' not in sys.path: sys.path.append('./')

from root.config import *
from SCREWS.frozen import FrozenOnly




class _3dCSCG_Edge(FrozenOnly):
    def __init__(self, mesh):
        self._mesh_ = mesh
        self._elements_ = None
        self._freeze_self_()


    @property
    def elements(self):
        """The edge elements. Only generate them when called first time."""
        if self._elements_ is None:
            self._elements_ = _3dCSCG_Edge_Elements(self)
        return self._elements_





class _3dCSCG_Edge_Elements(FrozenOnly):
    """"""
    def __init__(self, edge):
        """"""
        self._edge_ = edge
        self._mesh_ = edge._mesh_
        self._MAP_ = self.___PRIVATE_generating_edge_element_map___()
        self._locations_ = self.___PRIVATE_generating_edge_locations___()
        self._elements_ = dict()
        self._ct_ = _3dCSCG_Edge_Elements_CT(self)
        self._freeze_self_()

    def ___PRIVATE_generating_edge_element_map___(self):
        """"""
        MAP = dict()
        global_numbering = None
        # non-hybrid numbering ...
        mesh = self._mesh_
        if mesh.domain.IS_periodic:
            if rAnk == mAster_rank:
                baseElementLayout = mesh.elements.layout
                for rn in baseElementLayout:
                    regionElementLayout = baseElementLayout[rn]
                    assert all(np.array(regionElementLayout) > 1), \
                        f" elements.layout[{rn}]={regionElementLayout} wrong," \
                        f" needs (>1, >1, >1) to make it work for periodic domain."

        if rAnk != mAster_rank:
            element_map = mesh.elements.map
            element_indices = mesh.elements.indices
            cOmm.send([element_map, element_indices], dest=mAster_rank, tag=rAnk)
            global_numbering = cOmm.recv(source=mAster_rank, tag=rAnk)
        else:
            p = [1,1,1]
            numOfBasisComponents = [4,4,4]
            numberingCache = dict()
            currentNumber = 0
            other_side_name = 'SNEWFB' # not an error, this is other side name.
            sidePairDict = {'N':'S', 'S':'N', 'W':'E', 'E':'W', 'B':'F', 'F':'B'}
            for i in range(sIze):
                if i == mAster_rank:
                    element_map = mesh.elements.map
                    element_indices = mesh.elements.indices
                else:
                    element_map, element_indices = cOmm.recv(source=i, tag=i)

                for k in element_indices:
                    if k == 0:
                        assert numberingCache == dict(), "We make sure #0 element is numbered firstly."
                        nC0 = np.arange(currentNumber, currentNumber+numOfBasisComponents[0]).reshape(
                            (p[0], p[1]+1, p[2]+1), order='F')
                        currentNumber += numOfBasisComponents[0]
                        nC1 = np.arange(currentNumber, currentNumber+numOfBasisComponents[1]).reshape(
                            (p[0]+1, p[1], p[2]+1), order='F')
                        currentNumber += numOfBasisComponents[1]
                        nC2 = np.arange(currentNumber, currentNumber+numOfBasisComponents[2]).reshape(
                            (p[0]+1, p[1]+1, p[2]), order='F')
                        currentNumber += numOfBasisComponents[2]
                        numberingCache[0] = [nC0, nC1, nC2]
                        del nC0, nC1, nC2

                    else:
                        # TODO: see below explanation.
                        # note that this method is not good since it depends on the topology of the regions.
                        # When we number an edges of a mesh element #i, we need at least one of the six mesh
                        # elements that are face to face adjacent to mesh element #i to be already numbered.
                        # for some regions of special topology, this may not be the case, and it will raise
                        # exception since `k` in not in `numberingCache`. This also happens for the node
                        # elements and the Naive numbering of 0-form, 1-form, 2-form. But normally, this can
                        # be overcome with a proper sequence of regions indicated by `mesh.domain.regions.names`.
                        # we can define a specific sequence in domain_input initial function by setting
                        # `region_sequence`.
                        notNumberedPlaces0 = numberingCache[k][0] == -1
                        notNumberedPlaces1 = numberingCache[k][1] == -1
                        notNumberedPlaces2 = numberingCache[k][2] == -1
                        howManyNotNumbered0 = len(numberingCache[k][0][notNumberedPlaces0])
                        howManyNotNumbered1 = len(numberingCache[k][1][notNumberedPlaces1])
                        howManyNotNumbered2 = len(numberingCache[k][2][notNumberedPlaces2])
                        # assert howManyNotNumbered < numOfBasis, "We make sure at least a part is numbered!"
                        if howManyNotNumbered0 > 0:
                            numberingCache[k][0][notNumberedPlaces0] = np.arange(
                                currentNumber, currentNumber+howManyNotNumbered0
                            )
                            currentNumber += howManyNotNumbered0
                        if howManyNotNumbered1 > 0:
                            numberingCache[k][1][notNumberedPlaces1] = np.arange(
                                currentNumber, currentNumber+howManyNotNumbered1
                            )
                            currentNumber += howManyNotNumbered1
                        if howManyNotNumbered2 > 0:
                            numberingCache[k][2][notNumberedPlaces2] = np.arange(
                                currentNumber, currentNumber+howManyNotNumbered2
                            )
                            currentNumber += howManyNotNumbered2

                    for j, EMki in enumerate(element_map[k]):
                        if isinstance(EMki, str): # on domain boundary
                            pass
                        else:
                            otherElement = EMki
                            otherSide = other_side_name[j]
                            if otherElement not in numberingCache and otherElement > k:
                                numberingCache[otherElement] = [-np.ones((p[0], p[1]+1, p[2]+1), dtype=int),
                                                                -np.ones((p[0]+1, p[1], p[2]+1), dtype=int),
                                                                -np.ones((p[0]+1, p[1]+1, p[2]), dtype=int)]
                            if otherElement > k:
                                self.___PRIVATE_for_1Form_pass_element_side_numbering_from_to___(
                                    k, sidePairDict[otherSide], otherElement, otherSide, numberingCache
                                )

                toBeSentAway = dict()
                for k in element_indices: toBeSentAway[k] = numberingCache[k]
                if i == mAster_rank:
                    global_numbering = toBeSentAway
                else:
                    cOmm.send(toBeSentAway, dest=i, tag=i)
                for k in element_indices: del numberingCache[k]


        MAX = list()
        for i in self._mesh_.elements:
            vector = np.concatenate([global_numbering[i][0].ravel('F'),
                                     global_numbering[i][1].ravel('F'),
                                     global_numbering[i][2].ravel('F')])
            MAP[i] = list(vector)
            MAX.append(max(MAP[i]))

        if len(MAX) == 0:
            MAX = -1
        else:
            MAX = max(MAX)

        MAX = cOmm.gather(MAX, root=mAster_rank)
        if rAnk == mAster_rank:
            MAX = max(MAX) + 1
            assert MAX >= 4

        self._GLOBAL_num_ = cOmm.bcast(MAX, root=mAster_rank)

        return MAP

    @staticmethod
    def ___PRIVATE_for_1Form_pass_element_side_numbering_from_to___(
        fromElement, fromSide, toElement, toSide, numberingCache):

        data0, data1, data2 = None, None, None

        if fromSide == 'N'  :
            data1 = numberingCache[fromElement][1][0 , :, :]
            data2 = numberingCache[fromElement][2][0 , :, :]
        elif fromSide == 'S':
            data1 = numberingCache[fromElement][1][-1, :, :]
            data2 = numberingCache[fromElement][2][-1, :, :]
        elif fromSide == 'W':
            data0 = numberingCache[fromElement][0][ :, 0, :]
            data2 = numberingCache[fromElement][2][ :, 0, :]
        elif fromSide == 'E':
            data0 = numberingCache[fromElement][0][ :,-1, :]
            data2 = numberingCache[fromElement][2][ :,-1, :]
        elif fromSide == 'B':
            data0 = numberingCache[fromElement][0][ :, :, 0]
            data1 = numberingCache[fromElement][1][ :, :, 0]
        elif fromSide == 'F':
            data0 = numberingCache[fromElement][0][ :, :,-1]
            data1 = numberingCache[fromElement][1][ :, :,-1]
        else:
            raise Exception()

        if toSide == 'N'  :
            numberingCache[toElement][1][ 0, :, :] = data1
            numberingCache[toElement][2][ 0, :, :] = data2
        elif toSide == 'S':
            numberingCache[toElement][1][-1, :, :] = data1
            numberingCache[toElement][2][-1, :, :] = data2
        elif toSide == 'W':
            numberingCache[toElement][0][ :, 0, :] = data0
            numberingCache[toElement][2][ :, 0, :] = data2
        elif toSide == 'E':
            numberingCache[toElement][0][ :,-1, :] = data0
            numberingCache[toElement][2][ :,-1, :] = data2
        elif toSide == 'B':
            numberingCache[toElement][0][ :, :, 0] = data0
            numberingCache[toElement][1][ :, :, 0] = data1
        elif toSide == 'F':
            numberingCache[toElement][0][ :, :,-1] = data0
            numberingCache[toElement][1][ :, :,-1] = data1
        else:
            raise Exception()


    @property
    def GLOBAL_num(self):
        """How many edge elements in total (in all cores)?"""
        return self._GLOBAL_num_

    @property
    def map(self):
        """(dict) A dict whose keys are mesh element numbers and values
        are its 12 edge elements in the sequence ['WB', 'EB', 'WF', 'EF', 'NB', 'SB', 'NF', 'SF', 'NW', 'SW', 'NE', 'SE']."""
        return self._MAP_

    def ___PRIVATE_generating_edge_locations___(self):
        """"""
        MAP = self.map
        MAP = cOmm.gather(MAP, root=mAster_rank)
        if rAnk == mAster_rank:
            ___ = dict()
            for mp in MAP:
                add_len = len(mp)
                cur_len = len(___)
                ___.update(mp)
                assert len(___) == cur_len + add_len, "A trivial check."
            MAP = ___
            assert len(MAP) == self._mesh_.elements.GLOBAL_num, "A trivial check."

            LOC_DICT = dict()
            LOC_DICT_FULL = dict()
            ind_2_loc = ['WB', 'EB', 'WF', 'EF', 'NB', 'SB', 'NF', 'SF', 'NW', 'SW', 'NE', 'SE']

            for i in range(self._mesh_.elements.GLOBAL_num):
                mp_i = MAP[i]
                for ind, edge in enumerate(mp_i):
                    if edge not in LOC_DICT: LOC_DICT[edge] = list()
                    LOC_DICT[edge].append(str(i)+ind_2_loc[ind])

                    if len(LOC_DICT[edge]) == 4:
                        LOC_DICT_FULL[edge] = LOC_DICT[edge]
                        del LOC_DICT[edge]

            LOC_DICT_FULL.update(LOC_DICT) # LOC_DICT_FULL has all locations on mesh elements.
            del LOC_DICT

            # now we split the LOC_DICT_FULL to send to each cores.
            EDs = self._mesh_._element_distribution_
            LOCAL_EDGES = list()
            for core in EDs:
                local_elements = EDs[core]
                local_edges = dict()
                for i in local_elements:
                    mp_i = MAP[i]
                    for edge in mp_i:
                        if edge not in local_edges:
                            local_edges[edge] = LOC_DICT_FULL[edge]

                LOCAL_EDGES.append(local_edges)
            del LOC_DICT_FULL, local_edges

        else:
            LOCAL_EDGES = None

        LOCAL_EDGES = cOmm.scatter(LOCAL_EDGES, root=mAster_rank)
        mesh_map = self._mesh_.elements.map

        if saFe_mode:
            for i in self.map:
                for edge in self.map[i]:
                    assert edge in LOCAL_EDGES, f"safety check!"

        face_ind_dict = {'N':0, 'S':1, 'W':2, 'E':3, 'B':4, 'F':5}

        LOCAL_EDGES_BNS = dict()

        for edge in LOCAL_EDGES:
            for loc in LOCAL_EDGES[edge]:
                mesh_element = int(loc[:-2])
                corner = loc[-2:]
                if mesh_element in mesh_map: # we are looking at a location from a local mesh element.
                    for f in corner:
                        ind = face_ind_dict[f]
                        what_is_here = mesh_map[mesh_element][ind]
                        if isinstance(what_is_here, str): # a mesh boundary
                            if edge not in LOCAL_EDGES_BNS:
                                LOCAL_EDGES_BNS[edge] = list()
                            if what_is_here not in LOCAL_EDGES_BNS[edge]:
                                LOCAL_EDGES_BNS[edge].append(what_is_here)

        for edge in LOCAL_EDGES:
            if edge in LOCAL_EDGES_BNS:
                LOCAL_EDGES[edge].extend(LOCAL_EDGES_BNS[edge])

        for i in range(sIze):
            to_be_check = cOmm.bcast(LOCAL_EDGES, root=i)

            send_back = dict()
            if rAnk == i:
                pass
            else:
                for edge in to_be_check:

                    if edge in LOCAL_EDGES_BNS:
                        for loc in LOCAL_EDGES_BNS[edge]:
                            if loc not in to_be_check[edge]:
                                if edge not in send_back:
                                    send_back[edge] = list()
                                if loc not in send_back[edge]:
                                    send_back[edge].append(loc)

            send_back = cOmm.gather(send_back, root=i)

            if rAnk == i:
                for back in send_back:
                    for edge in back:
                        assert edge in LOCAL_EDGES, "safety check!"
                        for bn in back[edge]:
                            if bn not in LOCAL_EDGES[edge]:
                                LOCAL_EDGES[edge].append(bn)
            else:
                pass

        for i in self._MAP_:
            for edge in self._MAP_[i]:
                assert edge in LOCAL_EDGES, "safety check!"
                LOCAL_EDGES[edge] = tuple(LOCAL_EDGES[edge])

        bns = self._mesh_.boundaries.names
        ___ = '0123456789'

        shared_by_elements = dict()
        on_mesh_boundaries = dict()
        for edge in LOCAL_EDGES:
            shared_by_elements[edge] = list()
            on_mesh_boundaries[edge] = list()
            for position in LOCAL_EDGES[edge]:
                if position[0] in ___:
                    element = int(position[:-2])
                    if element not in shared_by_elements[edge]:
                        shared_by_elements[edge].append(element)
                else:
                    assert position in bns, f"position={position} is wrong."
                    assert position not in on_mesh_boundaries[edge], f"boundary: {position} appears more than once."
                    # noinspection PyUnresolvedReferences
                    on_mesh_boundaries[edge].append(position)

            # noinspection PyUnresolvedReferences
            shared_by_elements[edge] = tuple(shared_by_elements[edge]) # tuple is faster and memory less.
            # noinspection PyUnresolvedReferences
            on_mesh_boundaries[edge] = tuple(on_mesh_boundaries[edge]) # tuple is faster and memory less.

        self._shared_by_elements_ = shared_by_elements
        self._on_mesh_boundaries_ = on_mesh_boundaries

        PES = self._mesh_.___local_periodic_element_sides___
        LPE = self._mesh_.___local_periodic_elements___

        IS_on_periodic_boundary = dict()
        for edge in self._shared_by_elements_:
            IS_on_periodic_boundary[edge] = False
            if all([_ not in LPE for _ in self._shared_by_elements_[edge]]):
                pass
            else:
                for position in LOCAL_EDGES[edge]:
                    if position[0] in ___:
                        element, corners = position[:-2], position[-2:]
                        for corner in corners:
                            ele_side = element+corner
                            if ele_side in PES:
                                IS_on_periodic_boundary[edge] = True
                                break
                    else:
                        pass

                    if IS_on_periodic_boundary[edge]:
                        break

        self._IS_on_periodic_boundary_ = IS_on_periodic_boundary

        return LOCAL_EDGES

    @property
    def num(self):
        """Return how many local edge elements.
        (int)

        :return:
        """
        return len(self._locations_)

    def __getitem__(self, item):
        if item not in self._elements_:
            assert item in self._locations_, f"edge element #{item} is not included in this core (#{rAnk})."
            self._elements_[item] = _3dCSCG_Edge_Element(self, item)
        return self._elements_[item]

    def __contains__(self, item):
        return item in self._locations_

    def __iter__(self):
        for i in self._locations_:
            yield i

    def __len__(self):
        return self.num

    @staticmethod
    def ___generate_full_ep___(ep3, element_corner_edge):
        """

        :param ep3:
        :param element_corner_edge:
        :return:
        """
        assert len(ep3) == 3 and all([np.ndim(epi) == 1 for epi in ep3]), f"When we parse_1d_3ep, we need 3 evaluation_points of ndim=1 ."

        for i, ep in enumerate(ep3):
            assert np.max(ep) <= 1 and np.min(ep) >= -1 and np.all(
                np.diff(ep) > 0), f"evaluation_points[{i}] = {ep} is not an increasing 1d array in [-1,1]."

        ep0, ep1, ep2 = ep3
        s0 = np.shape(ep0)
        s1 = np.shape(ep1)
        s2 = np.shape(ep2)

        if element_corner_edge == 'WB':
            return ep0, -np.ones(s0), -np.ones(s0)
        elif element_corner_edge == 'EB':
            return ep0, np.ones(s0), -np.ones(s0)
        elif element_corner_edge == 'WF':
            return ep0, -np.ones(s0), np.ones(s0)
        elif element_corner_edge == 'EF':
            return ep0, np.ones(s0), np.ones(s0)
        elif element_corner_edge == 'NB':
            return -np.ones(s1), ep1, -np.ones(s1)
        elif element_corner_edge == 'SB':
            return np.ones(s1), ep1, -np.ones(s1)
        elif element_corner_edge == 'NF':
            return -np.ones(s1), ep1, np.ones(s1)
        elif element_corner_edge == 'SF':
            return np.ones(s1), ep1, np.ones(s1)
        elif element_corner_edge == 'NW':
            return -np.ones(s2), -np.ones(s2), ep2
        elif element_corner_edge == 'SW':
            return np.ones(s2), -np.ones(s2), ep2
        elif element_corner_edge == 'NE':
            return -np.ones(s2), np.ones(s2), ep2
        elif element_corner_edge == 'SE':
            return np.ones(s2), np.ones(s2), ep2
        else:
            raise Exception()

    @property
    def coordinate_transformation(self):
        """We will not cache the outputs because it is very fast anyway."""
        return self._ct_





class _3dCSCG_Edge_Elements_CT(FrozenOnly):
    """"""
    def __init__(self, elements):
        """"""
        self._elements_ = elements
        self._freeze_self_()











class _3dCSCG_Edge_Element(FrozenOnly):
    """"""
    def __init__(self, elements, i):
        """"""
        self._elements_ = elements
        self._i_ = i
        self._positions_ = elements._locations_[i]
        self._ct_ = None
        self._cp_ = None
        self._ce_ = None
        self._cce_ = None
        self._freeze_self_()

    @property
    def i(self):
        return self._i_

    @property
    def positions(self):
        """This edge element is at these positions."""
        return self._positions_

    @property
    def coordinate_transformation(self):
        if self._ct_ is None:
            self._ct_ = _3dCSCG_Edge_Element_CT(self)
        return self._ct_

    @property
    def shared_by_mesh_elements(self):
        return self._elements_._shared_by_elements_[self._i_]

    @property
    def on_mesh_boundaries(self):
        return self._elements_._on_mesh_boundaries_[self._i_]

    @property
    def IS_on_mesh_boundary(self):
        return True if len(self.on_mesh_boundaries) > 0 else False

    @property
    def IS_on_periodic_boundary(self):
        return self._elements_._IS_on_periodic_boundary_[self._i_]


    @property
    def CHARACTERISTIC_position(self):
        """The position we mainly locate this edge element."""
        if self._cp_ is None:

            for pos in self.positions:
                element, corner_edge = pos[:-2], pos[-2:] # before reach boundary names, we must have found it. So no worries.
                element = int(element)
                if element in self._elements_._MAP_:
                    self._cp_ = pos
                    self._ce_ = element
                    self._cce_ = corner_edge
                    break

        return self._cp_
    @property
    def CHARACTERISTIC_element(self):
        """We mainly consider this edge element is a corner-edge of this mesh
        element."""
        if self._ce_ is None:
            _ = self.CHARACTERISTIC_position
        return self._ce_
    @property
    def CHARACTERISTIC_corner_edge(self):
        """We main consider this edge element is such a corner-edge of the
        CHARACTERISTIC_element."""
        if self._cce_ is None:
            _ = self.CHARACTERISTIC_position
        return self._cce_
    @property
    def CHARACTERISTIC_region(self):
        """We mainly consider this edge element is in this region."""
        region = self._elements_._mesh_.DO.FIND_region_name_of_element(self.CHARACTERISTIC_element)
        return region






class _3dCSCG_Edge_Element_CT(FrozenOnly):
    """"""
    def __init__(self, ee):
        """"""
        self._ee_ = ee
        self._freeze_self_()

    def mapping(self, *ep3, from_element=None, corner_edge=None):
        """
        If `from_element` and `corner` are None, we compute it from this position.

        :param ep3:
        :param from_element: We compute it from this element.
        :param corner_edge: We compute it from this corner.
        :return:
        """
        if self._ee_.IS_on_periodic_boundary:
            assert from_element is not None, \
                "to compute the physical position of an edge element on periodic " \
                "boundary, we have to provide from which element you " \
                "want to compute it since it clearly will gives " \
                "different results."
            if from_element == 'any':
                # the different results do not matter; for example, when
                # we want to get value from a periodic function, the
                # location for evaluating the function also does not
                # matter.
                from_element = self._ee_.CHARACTERISTIC_element
            else:
                pass

        if from_element is None:
            i = self._ee_.CHARACTERISTIC_element
        elif from_element == 'any':
            i = self._ee_.CHARACTERISTIC_element
        else:
            i = from_element

        assert self._ee_.i in self._ee_._elements_.map[i], \
            f"edge element #{self._ee_.i} is not on mesh element {i}."

        ___ = ['WB', 'EB', 'WF', 'EF', 'NB', 'SB', 'NF', 'SF', 'NW', 'SW', 'NE', 'SE']
        if self._ee_._elements_.map[i].count(self._ee_.i) == 1: # this mesh element is not periodic to itself.
            corner_index = self._ee_._elements_.map[i].index(self._ee_.i)
            element_corner = ___[corner_index]
            if from_element is None: # if we do not provide `from_element` we must have this
                assert element_corner == self._ee_.CHARACTERISTIC_corner_edge

            if corner_edge is not None:
                assert element_corner == corner_edge, f"cannot compute it at provided corner edge: {corner_edge}"

        elif self._ee_._elements_.map[i].count(self._ee_.i) > 1: # this mesh element is periodic to itself.
            assert corner_edge is not None, f"edge element #{self._ee_.i} " \
                                            f"is on more than 1 corner-edges of element #{i} " \
                                            f"(periodic), provide corner_edge as well."
            element_corner = corner_edge
        else:
            raise Exception()

        # we will compute the physical position of this edge element from mesh element #`i` at its corner_edge `element_corner`
        assert self._ee_.i == self._ee_._elements_.map[i][___.index(element_corner)], \
            f"node element #{self._ee_.i} is not at {element_corner} of mesh element #{i}."

        ep = self._ee_._elements_.___generate_full_ep___(ep3, element_corner)
        x, y, z = self._ee_._elements_._mesh_.elements[i].coordinate_transformation.mapping(*ep)

        return x, y, z




if __name__ == '__main__':
    # mpiexec -n 12 python _3dCSCG\mesh\edge\main.py
    from _3dCSCG.main import MeshGenerator
    elements = [2, 2, 2]
    # mesh = MeshGenerator('crazy_periodic', c=0.0, bounds=([0,3], [0,3], [0,3]))(elements)
    mesh = MeshGenerator('bridge_arch_cracked')(elements)
    edges = mesh.edge.elements
