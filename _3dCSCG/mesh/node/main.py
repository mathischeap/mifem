# -*- coding: utf-8 -*-
"""
The edge elements of a mesh.
"""


import sys
if './' not in sys.path: sys.path.append('../')

from root.config import *
from screws.frozen import FrozenOnly





class _3dCSCG_Node(FrozenOnly):
    def __init__(self, mesh):
        self._mesh_ = mesh
        self._elements_ = None
        self._freeze_self_()


    @property
    def elements(self):
        """The edge elements. Only generate them when called first time."""
        if self._elements_ is None:
            self._elements_ = _3dCSCG_Node_Elements(self)
        return self._elements_





class _3dCSCG_Node_Elements(FrozenOnly):
    """"""
    def __init__(self, node):
        """"""
        self._node_ = node
        self._mesh_ = node._mesh_



        self._MAP_ = self.___PRIVATE_generating_node_element_map___()
        self._locations_ = self.___PRIVATE_generating_node_elements___()
        self._elements_ = dict()
        self._freeze_self_()

    def ___PRIVATE_generating_node_element_map___(self):
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
                        f" needs (>1, >1, >1) to make it work for periodic regions."



        if rAnk != mAster_rank:
            element_map = mesh.elements.map
            element_indices = mesh.elements.indices
            cOmm.send([element_map, element_indices], dest=mAster_rank, tag=rAnk)
            global_numbering = cOmm.recv(source=mAster_rank, tag=rAnk)
        else:
            p = (2, 2, 2)
            numberingCache = dict()
            numOfBasis = 8
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
                        numberingCache[0] = np.arange(numOfBasis).reshape(p, order='F')
                        currentNumber += numOfBasis
                    else:
                        notNumberedPlaces = numberingCache[k] == -1
                        howManyNotNumbered = len(numberingCache[k][notNumberedPlaces])
                        # assert howManyNotNumbered < numOfBasis, "We make sure at least a part is numbered!"
                        if howManyNotNumbered > 0:
                            numberingCache[k][notNumberedPlaces] = np.arange(
                                currentNumber, currentNumber+howManyNotNumbered
                            )
                            currentNumber += howManyNotNumbered

                    for j, EMki in enumerate(element_map[k]):
                        if isinstance(EMki, str): # on domain boundary
                            pass
                        else:
                            otherElement = EMki
                            otherSide = other_side_name[j]
                            if otherElement not in numberingCache and otherElement > k:
                                numberingCache[otherElement] = -np.ones(p, dtype=int)
                            if otherElement > k:
                                self.___PRIVATE_for_0Form_pass_element_side_numbering_from_to___(
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
            vector = global_numbering[i].ravel('F')
            MAP[i] = list(vector)
            MAX.append(max(MAP[i]))

        if len(MAX) == 0:
            MAX = -1
        else:
            MAX = max(MAX)

        MAX = cOmm.gather(MAX, root=mAster_rank)
        if rAnk == mAster_rank:
            MAX = max(MAX) + 1

        self._GLOBAL_num_ = cOmm.bcast(MAX, root=mAster_rank)

        return MAP

    @staticmethod
    def ___PRIVATE_for_0Form_pass_element_side_numbering_from_to___(
        fromElement, fromSide, toElement, toSide, numberingCache):
        if fromSide == 'N'  : data = numberingCache[fromElement][0 , :, :]
        elif fromSide == 'S': data = numberingCache[fromElement][-1, :, :]
        elif fromSide == 'W': data = numberingCache[fromElement][ :, 0, :]
        elif fromSide == 'E': data = numberingCache[fromElement][ :,-1, :]
        elif fromSide == 'B': data = numberingCache[fromElement][ :, :, 0]
        elif fromSide == 'F': data = numberingCache[fromElement][ :, :,-1]
        else: raise Exception()

        if toSide == 'N'  : numberingCache[toElement][ 0, :, :] = data
        elif toSide == 'S': numberingCache[toElement][-1, :, :] = data
        elif toSide == 'W': numberingCache[toElement][ :, 0, :] = data
        elif toSide == 'E': numberingCache[toElement][ :,-1, :] = data
        elif toSide == 'B': numberingCache[toElement][ :, :, 0] = data
        elif toSide == 'F': numberingCache[toElement][ :, :,-1] = data
        else: raise Exception()

    @property
    def map(self):
        """(dict) A dict whose keys are mesh element numbers and values
        are its 8 node elements in the sequence ['NWB', 'SWB', 'NEB', 'SEB', 'NWF', 'SWF', 'NEF', 'SEF']."""
        return self._MAP_

    def ___PRIVATE_generating_node_elements___(self):
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
            ind_2_loc = ['NWB', 'SWB', 'NEB', 'SEB', 'NWF', 'SWF', 'NEF', 'SEF']

            for i in range(self._mesh_.elements.GLOBAL_num):
                mp_i = MAP[i]
                for ind, node in enumerate(mp_i):
                    if node not in LOC_DICT: LOC_DICT[node] = list()
                    LOC_DICT[node].append(str(i)+ind_2_loc[ind])

                    if len(LOC_DICT[node]) == 8:
                        LOC_DICT_FULL[node] = LOC_DICT[node]
                        del LOC_DICT[node]

            LOC_DICT_FULL.update(LOC_DICT) # LOC_DICT_FULL has all locations on mesh elements.
            del LOC_DICT

            # now we split the LOC_DICT_FULL to send to each cores.
            EDs = self._mesh_._element_distribution_
            LOCAL_NODES = list()
            for core in EDs:
                local_elements = EDs[core]
                local_nodes = dict()
                for i in local_elements:
                    mp_i = MAP[i]
                    for node in mp_i:
                        if node not in local_nodes:
                            local_nodes[node] = LOC_DICT_FULL[node]

                LOCAL_NODES.append(local_nodes)
            del LOC_DICT_FULL, local_nodes

        else:
            LOCAL_NODES = None

        LOCAL_NODES = cOmm.scatter(LOCAL_NODES, root=mAster_rank)
        mesh_map = self._mesh_.elements.map

        if saFe_mode:
            for i in self.map:
                for node in self.map[i]:
                    assert node in LOCAL_NODES, f"safety check!"

        face_ind_dict = {'N':0, 'S':1, 'W':2, 'E':3, 'B':4, 'F':5}

        LOCAL_NODES_BNS = dict()

        for node in LOCAL_NODES:
            for loc in LOCAL_NODES[node]:
                mesh_element = int(loc[:-3])
                corner = loc[-3:]
                if mesh_element in mesh_map: # we are looking at a location from a local mesh element.
                    for f in corner:
                        ind = face_ind_dict[f]
                        what_is_here = mesh_map[mesh_element][ind]
                        if isinstance(what_is_here, str): # a mesh boundary
                            if node not in LOCAL_NODES_BNS:
                                LOCAL_NODES_BNS[node] = list()
                            if what_is_here not in LOCAL_NODES_BNS[node]:
                                LOCAL_NODES_BNS[node].append(what_is_here)

        for node in LOCAL_NODES:
            if node in LOCAL_NODES_BNS:
                LOCAL_NODES[node].extend(LOCAL_NODES_BNS[node])

        for i in range(sIze):
            to_be_check = cOmm.bcast(LOCAL_NODES, root=i)

            send_back = dict()
            if rAnk == i:
                pass
            else:
                for node in to_be_check:

                    if node in LOCAL_NODES_BNS:
                        for loc in LOCAL_NODES_BNS[node]:
                            if loc not in to_be_check[node]:
                                if node not in send_back:
                                    send_back[node] = list()
                                if loc not in send_back[node]:
                                    send_back[node].append(loc)


            send_back = cOmm.gather(send_back, root=i)

            if rAnk == i:
                for back in send_back:
                    for node in back:
                        assert node in LOCAL_NODES, "safety check!"
                        for bn in back[node]:
                            if bn not in LOCAL_NODES[node]:
                                LOCAL_NODES[node].append(bn)

            else:
                pass

        for i in self._MAP_:
            for node in self._MAP_[i]:
                assert node in LOCAL_NODES, "safety check!"
                LOCAL_NODES[node] = tuple(LOCAL_NODES[node])



        bns = self._mesh_.boundaries.names
        ___ = '0123456789'

        shared_by_elements = dict()
        on_mesh_boundaries = dict()
        for node in LOCAL_NODES:
            shared_by_elements[node] = list()
            on_mesh_boundaries[node] = list()
            for position in LOCAL_NODES[node]:
                if position[0] in ___:
                    element = int(position[:-3])
                    if element not in shared_by_elements[node]:
                        shared_by_elements[node].append(element)
                else:
                    assert position in bns, f"position={position} is wrong."
                    assert position not in on_mesh_boundaries[node], f"boundary: {position} appears more than once."
                    # noinspection PyUnresolvedReferences
                    on_mesh_boundaries[node].append(position)

            # noinspection PyUnresolvedReferences
            shared_by_elements[node] = tuple(shared_by_elements[node]) # tuple is faster and memory less.
            # noinspection PyUnresolvedReferences
            on_mesh_boundaries[node] = tuple(on_mesh_boundaries[node]) # tuple is faster and memory less.

        self._shared_by_elements_ = shared_by_elements
        self._on_mesh_boundaries_ = on_mesh_boundaries

        PES = self._mesh_.___local_periodic_element_sides___
        LPE = self._mesh_.___local_periodic_elements___

        IS_on_periodic_boundary = dict()
        for node in self._shared_by_elements_:
            IS_on_periodic_boundary[node] = False
            if all([_ not in LPE for _ in self._shared_by_elements_[node]]):
                pass
            else:
                for position in LOCAL_NODES[node]:
                    if position[0] in ___:
                        element, corners = position[:-3], position[-3:]
                        for corner in corners:
                            ele_side = element+corner
                            if ele_side in PES:
                                IS_on_periodic_boundary[node] = True
                                break
                    else:
                        pass

                    if IS_on_periodic_boundary[node]:
                        break
        self._IS_on_periodic_boundary_ = IS_on_periodic_boundary



        return LOCAL_NODES

    @property
    def GLOBAL_num(self):
        """How many node elements in total (in all cores)?"""
        return self._GLOBAL_num_

    @property
    def num(self):
        """Return how many local trace elements.
        (int)

        :return:
        """
        return len(self._locations_)

    def __getitem__(self, item):
        if item not in self._elements_:
            assert item in self._locations_, f"node element #{item} is not included in this core (#{rAnk})."
            self._elements_[item] = _3dCSCG_Node_Element(self, item)
        return self._elements_[item]

    def __contains__(self, item):
        return item in self._locations_

    def __iter__(self):
        for i in self._locations_:
            yield i

    def __len__(self):
        return self.num






class _3dCSCG_Node_Element(FrozenOnly):
    """"""
    def __init__(self, elements, i):
        """"""
        self._elements_ = elements
        self._i_ = i
        self._positions_ = elements._locations_[i]
        self._ct_ = None
        self._cp_ = None
        self._ce_ = None
        self._cc_ = None
        self._freeze_self_()

    @property
    def i(self):
        return self._i_

    @property
    def positions(self):
        """This node element is at these positions."""
        return self._positions_

    @property
    def coordinate_transformation(self):
        if self._ct_ is None:
            self._ct_ = _3dCSCG_Node_Element_CT(self)
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
        """The position we mainly locate this node element."""
        if self._cp_ is None:

            for pos in self.positions:
                element, corner = pos[:-3], pos[-3:] # before reach boundary names, we must have found it. So no worries.
                element = int(element)
                if element in self._elements_._MAP_:
                    self._cp_ = pos
                    self._ce_ = element
                    self._cc_ = corner
                    break
        return self._cp_
    @property
    def CHARACTERISTIC_element(self):
        """We mainly consider this node element is a corner of this mesh
        element."""
        if self._ce_ is None:
            _ = self.CHARACTERISTIC_position
        return self._ce_
    @property
    def CHARACTERISTIC_corner(self):
        """We mainly consider this node element is such a corner of the
        CHARACTERISTIC_element."""
        if self._cc_ is None:
            _ = self.CHARACTERISTIC_position
        return self._cc_
    @property
    def CHARACTERISTIC_region(self):
        """We mainly consider this node element is in this regions."""
        region = self._elements_._mesh_.DO.FIND_region_name_of_element(self.CHARACTERISTIC_element)
        return region



class _3dCSCG_Node_Element_CT(FrozenOnly):
    """"""
    def __init__(self, ne):
        """"""
        self._ne_ = ne
        self._freeze_self_()

    def mapping(self, from_element=None, corner=None):
        """
        If `from_element` and `corner` are None, we compute it from this position.

        :param from_element: We compute it from this element.
        :param corner: We compute it from this corner.
        :return:
        """
        if self._ne_.IS_on_periodic_boundary:
            assert from_element is not None, \
                "to compute the physical position of a node element on periodic " \
                "boundary, we have to provide from which element you " \
                "want to compute it since it clearly will gives " \
                "different results."
            if from_element == 'any':
                # the different results do not matter; for example, when
                # we want to get value from a periodic function, the
                # location for evaluating the function also does not
                # matter.
                from_element = self._ne_.CHARACTERISTIC_element
            else:
                pass


        if from_element is None:
            i = self._ne_.CHARACTERISTIC_element
        elif from_element == 'any':
            i = self._ne_.CHARACTERISTIC_element
        else:
            i = from_element


        assert self._ne_.i in self._ne_._elements_.map[i], \
            f"trace element #{self._ne_.i} is not on mesh element {i}."

        ___ = ['NWB', 'SWB', 'NEB', 'SEB', 'NWF', 'SWF', 'NEF', 'SEF']
        if self._ne_._elements_.map[i].count(self._ne_.i) == 1: # this mesh element is not periodic to itself.
            corner_index = self._ne_._elements_.map[i].index(self._ne_.i)
            element_corner = ___[corner_index]
            if from_element is None: # if we do not provide `from_element` we must have this
                assert element_corner == self._ne_.CHARACTERISTIC_corner

            if corner is not None:
                assert element_corner == corner, f"cannot compute it at provided corner {corner}"

        elif self._ne_._elements_.map[i].count(self._ne_.i) > 1: # this mesh element is periodic to itself.
            assert corner is not None, f"node element #{self._ne_.i} " \
                                       f"is on more than 1 corners of element #{i} " \
                                       f"(periodic), provide corner as well."
            element_corner = corner
        else:
            raise Exception()

        # we will compute the physical position of this node element from mesh element #`i` at its corner `element_corner`
        assert self._ne_.i == self._ne_._elements_.map[i][___.index(element_corner)], \
            f"node element #{self._ne_.i} is not at {element_corner} of mesh element #{i}."

        ep = self.___generate_full_ep___(element_corner)
        x, y, z = self._ne_._elements_._mesh_.elements[i].coordinate_transformation.mapping(*ep)

        return x, y, z

    @staticmethod
    def ___generate_full_ep___(element_corner):
        """
        :param element_corner:
        :return:
        """
        if element_corner == 'NWB':
            return -1, -1, -1
        elif element_corner == 'SWB':
            return 1, -1, -1
        elif element_corner == 'NEB':
            return -1, 1, -1
        elif element_corner == 'SEB':
            return 1, 1, -1
        elif element_corner == 'NWF':
            return -1, -1, 1
        elif element_corner == 'SWF':
            return 1, -1, 1
        elif element_corner == 'NEF':
            return -1, 1, 1
        elif element_corner == 'SEF':
            return 1, 1, 1
        else:
            raise Exception()














if __name__ == '__main__':
    # mpiexec -n 12 python _3dCSCG\mesh\node.py
    from _3dCSCG.main import MeshGenerator
    elements = [2, 2, 2]
    mesh = MeshGenerator('crazy_periodic', c=0.0, bounds=([0,3], [0,3], [0,3]))(elements)
    # mesh = MeshGenerator('bridge_arch_cracked')(elements)
    nodes = mesh.node.elements

    # print(rAnk, mesh.___local_periodic_element_sides___)

    print(nodes.GLOBAL_num)

    # for i in nodes:
    #     node = nodes[i]
        # print(i, edge.i, edge.positions, edge.CHARACTERISTIC_element, edge.CHARACTERISTIC_position)
        # print(i, node.i, node.positions, node.CHARACTERISTIC_position, node.CHARACTERISTIC_element, node.CHARACTERISTIC_corner, node.CHARACTERISTIC_region)
