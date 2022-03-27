


from root.config.main import *

from screws.freeze.main import FrozenOnly
from objects.CSCG._2d.mesh.trace.elements.element.main import _2dCSCG_Trace_Element



class _2dCSCG_Trace_Elements(FrozenOnly):
    def __init__(self, trace):
        self._trace_ = trace
        self._mesh_ = trace._mesh_
        self.___PRIVATE_generating_trace_elements___()
        self.___PRIVATE_reset_cache___()
        self._freeze_self_()


    def ___PRIVATE_reset_cache___(self):
        self._type_amount_dict_ = None

    def ___PRIVATE_find_type_and_amount_numbered_before___(self):
        """
        :return: A dictionary. For example, ``{..., 65: [32, 33,], ...}``, it means
            we have 32 'UD', 33 'LR' trace elements numbered before 65. We can see that
            32 + 33 = 65 (0, 1, 2, ..., 64).
        :rtype: dict
        """
        if self._type_amount_dict_ is None:
            if rAnk == 0:
                POOL = dict()
            else:
                POOL = cOmm.recv(source=rAnk - 1, tag=rAnk)

            type_amount_dict = dict()

            INDICES = list(self._elements_.keys())
            INDICES = sorted(INDICES)

            for i in INDICES:
                if i == 0: POOL[0] = np.array([0,0])

                tei = self[i]
                cs = tei.CHARACTERISTIC_edge
                POOL_i_p_1 = POOL[i].copy()
                if cs in 'UD':
                    POOL_i_p_1[0] += 1
                elif cs in 'LR':
                    POOL_i_p_1[1] += 1
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

                if not tei.IS_shared_by_cores:
                    del POOL[i]
                else:
                    CORE = tei.shared_with_core
                    if CORE < rAnk:
                        del POOL[i]

            if rAnk == sIze - 1:
                pass
            else:
                cOmm.send(POOL, dest=rAnk+1, tag=rAnk+1)

            self._type_amount_dict_ = type_amount_dict

        return self._type_amount_dict_


    @property
    def map(self):
        return self._MAP_

    @property
    def num(self):
        """
        (int)

        :return:
        """
        return len(self._elements_)

    def __getitem__(self, item):
        return self._elements_[item]

    def __contains__(self, item):
        return item in self._elements_

    def __iter__(self):
        for i in self._elements_:
            yield i

    def __len__(self):
        return self.num


    @property
    def GLOBAL_num(self):
        """
        (int) The total number of trace elements.

        :return:
        """
        return self._GLOBAL_num_

    def ___PRIVATE_generating_trace_elements___(self):
        self._elements_ = dict()
        elements_map = self._mesh_.elements.map
        sideNames = 'UDLR'
        sidePairs = {'U':'D', 'D':'U', 'L':'R', 'R':'L'}
        upesp = self._mesh_.___useful_periodic_element_edge_pairs___
        self._MAP_ = dict()

        if rAnk == 0: # first core
            cn = 0
            POOL = dict()
        else: # intermediate cores
            POOL, cn = cOmm.recv(source=rAnk-1, tag=rAnk)

        LOCAL_POOL = dict()
        for i in elements_map:
            self._MAP_[i] = [None for _ in range(4)]
            for j in range(4):
                side_1 = sideNames[j]
                position_1 = str(i) + side_1
                what = elements_map[i][j]
                if isinstance(what, str): #on domain boundary
                    position_2 = what
                    self._elements_[cn] = _2dCSCG_Trace_Element(
                        self, cn, position_1, position_2, position_1, ondb=True, onpb=False)
                    self._MAP_[i][j] = cn
                    cn += 1
                else:
                    side_2 = sidePairs[sideNames[j]]
                    position_2 = str(what) + side_2
                    pair_1 = str(i) + '-' + side_1 + '|' + side_2 + '-' + str(what)
                    pair_2 = str(what) + '-' + side_2 + '|' + side_1 + '-' + str(i)
                    onpb = True if pair_1 in upesp or pair_2 in upesp else False
                    if what > i:
                        self._elements_[cn] = _2dCSCG_Trace_Element(
                            self, cn, position_1, position_2, position_1, ondb=False, onpb=onpb)
                        self._MAP_[i][j] = cn
                        POOL[position_1] = cn
                        POOL[position_2] = cn
                        cn += 1
                    elif what == i:
                        assert onpb
                        if side_1 in 'UL':
                            self._elements_[cn] = _2dCSCG_Trace_Element(
                                self, cn, position_1, position_2, position_1, ondb=False, onpb=True)
                            self._MAP_[i][j] = cn
                            LOCAL_POOL[position_2] = cn
                            cn += 1
                        else:
                            self._MAP_[i][j] = LOCAL_POOL[position_1]
                            del LOCAL_POOL[position_1]
                    else:
                        assert position_1 in POOL and position_2 in POOL
                        assert POOL[position_1] == POOL[position_2]
                        tn = POOL[position_1]
                        self._MAP_[i][j] = tn
                        if tn in self._elements_:
                            pass
                        else:
                            self._elements_[tn] = _2dCSCG_Trace_Element(
                                self, tn, position_2, position_1, position_1, ondb=False, onpb=onpb)

                        del POOL[position_1]
                        del POOL[position_2]

        if rAnk == sIze - 1:
            assert len(POOL) == 0
        else:
            cOmm.send([POOL, cn], dest=rAnk+1, tag=rAnk+1)

        cOmm.barrier()

        cn = cOmm.bcast(cn, root=sIze-1)
        self._GLOBAL_num_ = cn

        if saFe_mode:
            for i in elements_map:
                for j in range(4):
                    assert self._MAP_[i][j] in self, f"Miss local trace element {self._MAP_[i][j]}"


    def DO_compute_mapping_of_trace_at_position(self, position, ep1d):
        """

        :param position: Like "10-U", then we compute for the trace element which is the upper edge of
            10th mesh element.
        :param ep1d: Just like the ep1d in element.ct.
        :return: Only return to master core.
        """
        edge = position[-1]
        element = int(position[:-1])

        if element in self._mesh_.elements:
            ep = self.___generate_full_ep___(ep1d, edge)
            xy = self._mesh_.elements[element].coordinate_transformation.mapping(*ep)
        else:
            xy = None

        xy = cOmm.gather(xy, root=mAster_rank)
        if rAnk == mAster_rank:
            assert xy.count(None) == sIze - 1
            return next(item for item in xy if item is not None) # find the first item which is not None.

    @staticmethod
    def ___generate_full_ep___(ep1d, element_edge):
        """
        Given the 1d evaluation points, we generate the 2d coordinates according to
        the element_edge.

        Parameters
        ----------
        ep1d:
            The 1d evaluation points.
        element_edge:
            The mesh-element edge: U, D, L, or R.

        """
        assert np.ndim(ep1d) == 1, " <TraceElementCoordinateTransformation> "
        assert np.min(ep1d) >= -1 and np.max(ep1d) <= 1, " <TraceElementCoordinateTransformation> "
        if np.size(ep1d) > 1:
            assert all(np.diff(ep1d) > 0), " <TraceElementCoordinateTransformation> "
        shape = len(ep1d)
        if element_edge == 'U':
            ep = (-np.ones(shape), ep1d)
        elif element_edge == 'D':
            ep = (np.ones(shape), ep1d)
        elif element_edge == 'L':
            ep = (ep1d, -np.ones(shape))
        elif element_edge == 'R':
            ep = (ep1d, np.ones(shape))
        else:
            raise Exception()
        return ep
