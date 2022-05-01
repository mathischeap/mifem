



import sys
if './' not in sys.path: sys.path.append('./')

from root.config.main import *
from screws.freeze.main import FrozenOnly


from objects.CSCG._3d.mesh.trace.elements.do.main import _3dCSCG_Trace_Elements_DO
from objects.CSCG._3d.mesh.trace.elements.selfcheck import _3dCSCG_Trace_Elements_SELFCHECK
from objects.CSCG._3d.mesh.trace.elements.group import _3dCSCG_Trace_Elements_Group
from objects.CSCG._3d.mesh.trace.elements.coordinate_transformation.main import _3dCSCG_Trace_Elements_CoordinateTransformation
from objects.CSCG._3d.mesh.trace.elements.element.main import _3dCSCG_Trace_Element



class _3dCSCG_Trace_Elements(FrozenOnly):
    """Note that a trace element can be in two cores because we store
    all 6 sides of each local mesh element. So if a trace element is on
    the interface of two mesh elements being in different cores, this
    trace element will be stored in two cores."""
    def __init__(self, trace):
        self._trace_ = trace
        self._mesh_ = trace._mesh_
        self._group_ = None
        self._DO_ = None
        self._SELFCHECK_ = None
        self._ct_ = None
        self.___NS___ = 'NS'
        self.___WE___ = 'WE'
        self.___BF___ = 'BF'
        self.___PRIVATE_generating_trace_elements___()
        self.___PRIVATE_reset_cache___()
        self._freeze_self_()


    def ___PRIVATE_reset_cache___(self):
        self._type_amount_dict_ = None
        self._Q_ = None
        self._AvQ_ = None
        self._WorstQ_ = None
        self._BestQ_ = None
        self.___multi_elements_metric___ = None


    def ___PRIVATE_find_type_and_amount_numbered_before___(self):
        """
        :return: A dictionary. For example, ``{..., 107: [32, 33, 42], ...}``, it means
            we have 32 'NS', 33 'WE' and 42 'BF' trace elements numbered before 107. We can see that
            32 + 33 + 42 = 107 (0, 1, 2, ..., 106).
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
                if i == 0: POOL[0] = np.array([0,0,0])

                tei = self[i]
                cs = tei.CHARACTERISTIC_side
                POOL_i_p_1 = POOL[i].copy()
                if cs in 'NS':
                    POOL_i_p_1[0] += 1
                elif cs in 'WE':
                    POOL_i_p_1[1] += 1
                elif cs in 'BF':
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

                if not tei.IS.shared_by_cores:
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
        """(dict) A dict whose keys are mesh element numbers and values
        are its 6 trace elements in the sequence N->S->W->E->B->F."""
        return self._MAP_

    @property
    def num(self):
        """Return how many local trace elements.
        (int)

        :return:
        """
        return len(self._elements_)

    @property
    def group(self):
        if self._group_ is None:
            self._group_ = _3dCSCG_Trace_Elements_Group(self)
        return self._group_

    @property
    def do(self):
        if self._DO_ is None:
            self._DO_ = _3dCSCG_Trace_Elements_DO(self)
        return self._DO_

    @property
    def selfcheck(self):
        if self._SELFCHECK_ is None:
            self._SELFCHECK_ = _3dCSCG_Trace_Elements_SELFCHECK(self)
        return self._SELFCHECK_

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
    def quality(self):
        """

        :return: A dict whose keys are the numbers of local trace
            elements and values are the qualities of the trace elements.

            quality will be [0,1]: 1, the best, 0: the worst.
        """
        if self._Q_ is None:
            self._Q_, self._AvQ_, self._WorstQ_, self._BestQ_ = \
                self.do.get_quality_of_trace_elements()
        return self._Q_

    @property
    def coordinate_transformation(self):
        if self._ct_ is None:
            self._ct_ = _3dCSCG_Trace_Elements_CoordinateTransformation(self)
        return self._ct_

    @property
    def _multi_elements_metric_(self):
        """Return a dict whose keys are the mark of type_wrt_metric of trace elements. And the values
        are the amount of trace elements that possess this type_wrt_metric."""
        if self.___multi_elements_metric___ is None:
            self.___multi_elements_metric___ = dict()

            for i in self:
                te = self[i]
                twm = te.type_wrt_metric.mark
                if isinstance(twm, int): # chaotic trace element, ignore it.
                    pass
                else:
                    if twm in self.___multi_elements_metric___:
                        self.___multi_elements_metric___[twm] += 1
                    else:
                        self.___multi_elements_metric___[twm] = 1

        return self.___multi_elements_metric___

    @property
    def GLOBAL_num(self):
        """
        (int) The total number of trace elements across all cores.

        :return:
        """
        return self._GLOBAL_num_

    def ___PRIVATE_generating_trace_elements___(self):
        self._elements_ = dict()
        elements_map = self._mesh_.elements.map
        sideNames = 'NSWEBF'
        sidePairs = {'N':'S', 'S':'N', 'W':'E', 'E':'W', 'B':'F', 'F':'B'}
        upesp = self._mesh_.___useful_periodic_element_side_pairs___
        self._MAP_ = dict()

        if rAnk == 0: # first core
            cn = 0
            POOL = dict()
        else: # intermediate cores
            POOL, cn = cOmm.recv(source=rAnk-1, tag=rAnk)

        LOCAL_POOL = dict()
        for i in elements_map:
            self._MAP_[i] = [None for _ in range(6)]
            for j in range(6):
                side_1 = sideNames[j]
                position_1 = str(i) + side_1
                what = elements_map[i][j]
                if isinstance(what, str): #on domain boundary
                    position_2 = what
                    self._elements_[cn] = _3dCSCG_Trace_Element(
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
                        self._elements_[cn] = _3dCSCG_Trace_Element(
                            self, cn, position_1, position_2, position_1, ondb=False, onpb=onpb)
                        self._MAP_[i][j] = cn
                        POOL[position_1] = cn
                        POOL[position_2] = cn
                        cn += 1
                    elif what == i:
                        assert onpb
                        if side_1 in 'NWB':
                            self._elements_[cn] = _3dCSCG_Trace_Element(
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
                            self._elements_[tn] = _3dCSCG_Trace_Element(
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
            map_set = set()
            for i in elements_map:
                map_set.update(self._MAP_[i])
            assert map_set == set(self._elements_.keys())
            assert self._MAP_.keys() == elements_map.keys()
            # should only contain trace elements which are sides of local mesh element.

    @staticmethod
    def ___generate_full_ep___(evaluation_points, element_side, parse_3_1d_eps=False, picking=False):
        """

        :param evaluation_points:
        :param element_side:
        :param parse_3_1d_eps: If `parse_ep` is True, then we have *ep is xi, eta, sigma, and they are all 1d,
            between [-1,1], we will pick up two from them according the side and do the mesh grid.

            When `parse_3_1d_eps`, we will automatically do the picking! So `picking` must be False.

        :param picking: We have provided 3 eps, but we will select two of them and make the third one
            according to the `element_side`.
        """
        if parse_3_1d_eps: # we actually received 3 1d arrays, and we will do meshgrid to make it proper.

            assert not picking, f"When `parse_3_1d_eps`, `picking` must be False " \
                                f"because we will do picking automatically."

            assert len(evaluation_points) == 3 and all([np.ndim(epi)==1 for epi in evaluation_points]), \
                f"When we parse_1d_3ep, we need 3 evaluation_points of ndim=1 ."

            for i, ep in enumerate(evaluation_points):
                assert np.max(ep) <= 1 and np.min(ep) >= -1 and np.all(np.diff(ep)>0), \
                    f"evaluation_points[{i}] = {ep} is not an increasing 1d array in [-1,1]."

            # evaluation_points will be of length 2.
            if element_side in 'NS':
                evaluation_points = np.meshgrid(evaluation_points[1], evaluation_points[2], indexing='ij')
            elif element_side in 'WE':
                evaluation_points = np.meshgrid(evaluation_points[0], evaluation_points[2], indexing='ij')
            elif element_side in 'BF':
                evaluation_points = np.meshgrid(evaluation_points[0], evaluation_points[1], indexing='ij')
            else:
                raise Exception()

        else: # we have standard evaluation_points, so do nothing here.
            pass

        if len(evaluation_points) == 2: # only two valid ep for local trace element provided.

            assert not picking, f"We already make the picking basically!"

            ep0, ep1 = evaluation_points
            shape0 = np.shape(ep0)

            if saFe_mode:
                assert shape0 == np.shape(ep1), \
                    " <TraceElement3D> : evaluation_points shape wrong."

            if   element_side == 'N': ep = (-np.ones(shape0), *evaluation_points)
            elif element_side == 'S': ep = ( np.ones(shape0), *evaluation_points)
            elif element_side == 'W': ep = (ep0, -np.ones(shape0), ep1)
            elif element_side == 'E': ep = (ep0,  np.ones(shape0), ep1)
            elif element_side == 'B': ep = (*evaluation_points, -np.ones(shape0))
            elif element_side == 'F': ep = (*evaluation_points,  np.ones(shape0))
            else:
                raise Exception()

            return ep

        elif len(evaluation_points) == 3:

            if picking: # we need to select from the 3 evaluation_points according to the element side.

                ep0, ep1, ep2 = evaluation_points
                shape0 = np.shape(ep0)
                shape1 = np.shape(ep1)
                shape2 = np.shape(ep2)
                assert shape0 == shape1 == shape2, \
                    f"When we do picking, for safety reason, " \
                    f"we ask all eps to be of same shape."

                if   element_side == 'N': ep = (-np.ones(shape0), ep1, ep2)
                elif element_side == 'S': ep = ( np.ones(shape0), ep1, ep2)
                elif element_side == 'W': ep = (ep0, -np.ones(shape1), ep2)
                elif element_side == 'E': ep = (ep0,  np.ones(shape1), ep2)
                elif element_side == 'B': ep = (ep0, ep1, -np.ones(shape2))
                elif element_side == 'F': ep = (ep0, ep1,  np.ones(shape2))
                else: raise Exception()

                return ep

            else: # two valid plus one -1 or +1 coordinates are provided
                if saFe_mode:
                    # noinspection PyUnresolvedReferences
                    assert evaluation_points[0].shape == \
                           evaluation_points[1].shape == \
                           evaluation_points[2].shape, "evaluation_points shape wrong."
                return evaluation_points

        else:
            raise Exception("evaluation_points shape wrong dimension wrong.")












if __name__ == '__main__':
    # mpiexec -n 12 python _3dCSCG\mesh\trace\elements\main.py
    from objects.CSCG._3d.master import MeshGenerator
    elements = [3, 4, 2]
    mesh = MeshGenerator('crazy_periodic', c=0.3, bounds=([0,1], [0,1], [0,1]))(elements)
    mesh.trace.elements.selfcheck.outward_unit_normal_vector()
    Q = mesh.trace.elements.quality
    # print(mesh.quality)
    # print(mesh.trace.quality)

    # mesh.trace.elements.do.illustrate_element(1)

    # te0 = mesh.trace.elements[0]

    # print(te0.IS.on_periodic_boundary)

    # for i in range(mesh.trace.elements.GLOBAL_num):
    #     mesh.trace.elements.do.illustrate_element(i)
        # if i in mesh.trace.elements:
        #     te = mesh.trace.elements[i]
        #
        #     print(rAnk, te.type_wrt_metric.mark)