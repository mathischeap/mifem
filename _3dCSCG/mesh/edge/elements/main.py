


import sys
if './' not in sys.path: sys.path.append('./')




from SCREWS.frozen import FrozenOnly
from root.config import *


from _3dCSCG.mesh.edge.elements.element.main import _3dCSCG_Edge_Element
from _3dCSCG.mesh.edge.elements.coordinate_transformation import _3dCSCG_Edge_Elements_CT


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
        self.RESET_cache()
        self._freeze_self_()

    def RESET_cache(self):
        pass

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
                                self.___PRIVATE_EDGE_for_1Form_pass_element_side_numbering_from_to___(
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
    def ___PRIVATE_EDGE_for_1Form_pass_element_side_numbering_from_to___(
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


    def ___DO_find_type_and_amount_numbered_before___(self):
        """
        :return: A dictionary. For example, ``{..., 107: [32, 33, 42], ...}``, it means
            we have 32 'NS'-, 33 'WE'- and 42 'BF'-direction edge elements numbered before 107.
            We can see that 32 + 33 + 42 = 107 (0, 1, 2, ..., 106).
        :rtype: dict
        """
        if rAnk == 0:
            POOL = dict()
        else:
            POOL = cOmm.recv(source=rAnk - 1, tag=rAnk)

        type_amount_dict = dict()

        INDICES = list(self._locations_.keys())
        INDICES = sorted(INDICES)

        for i in INDICES:
            if i == 0: POOL[0] = np.array([0,0,0])

            ee = self[i]
            cs = ee.direction
            POOL_i_p_1 = POOL[i].copy()
            if cs == 'NS':
                POOL_i_p_1[0] += 1
            elif cs == 'WE':
                POOL_i_p_1[1] += 1
            elif cs == 'BF':
                POOL_i_p_1[2] += 1
            else:
                raise Exception()
            if i+1 in POOL:
                assert np.all(POOL_i_p_1 == POOL[i+1])
            else:
                if i == self.GLOBAL_num - 1:
                    assert np.sum(POOL_i_p_1) == self.GLOBAL_num, "Something is wrong."
                else:
                    POOL[i+1] = POOL_i_p_1

            type_amount_dict[i] = POOL[i]

        if rAnk == sIze - 1:
            pass
        else:
            cOmm.send(POOL, dest=rAnk+1, tag=rAnk+1)

        for i in type_amount_dict:
            assert np.sum(type_amount_dict[i]) == i
            if i - 1 in type_amount_dict:
                A, B, C = type_amount_dict[i]
                a, b, c = type_amount_dict[i-1]
                if self[i-1].direction == 'NS':
                    assert A == a + 1 and B == b and C == c
                elif self[i-1].direction == 'WE':
                    assert A == a and B == b + 1 and C == c
                elif self[i-1].direction == 'BF':
                    assert A == a and B == b and C == c + 1
                else:
                    raise Exception()

        return type_amount_dict




    @property
    def GLOBAL_num(self):
        """How many edge elements in total (in all cores)?"""
        return self._GLOBAL_num_

    @property
    def map(self):
        """(dict) A dict whose keys are mesh element numbers and values
        are its 12 edge elements in the sequence ['WB', 'EB', 'WF', 'EF', 'NB', 'SB',
        'NF', 'SF', 'NW', 'SW', 'NE', 'SE'].
        """
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







if __name__ == '__main__':
    # mpiexec -n 12 python _3dCSCG\mesh\edge\elements\element\main.py
    from _3dCSCG.main import MeshGenerator
    elements = [2, 2, 2]
    # mesh = MeshGenerator('crazy_periodic', c=0.0, bounds=([0,3], [0,3], [0,3]))(elements)
    mesh = MeshGenerator('bridge_arch_cracked')(elements)
    edges = mesh.edge.elements

    for i in edges:
        edge = edges[i]