# -*- coding: utf-8 -*-
"""
The trace (face) elements of a mesh.
"""

import sys
if './' not in sys.path: sys.path.append('./')

from root.config import *
from SCREWS.frozen import FrozenOnly
from SCREWS.customized_warnings import TraceElementWarning
from SCREWS.functions._3d import angle_between_two_vectors
from itertools import combinations

import warnings
import matplotlib.pyplot as plt


class _3dCSCG_Trace(FrozenOnly):
    def __init__(self, mesh):
        self._mesh_ = mesh
        self._elements_ = _3dCSCG_Trace_Elements(self) # please initialize it here!
        self.RESET_cache()
        self._freeze_self_()

    def RESET_cache(self):
        pass

    @property
    def elements(self):
        return self._elements_

    @property
    def quality(self):
        _ = self.elements.quality
        AvQ, WorstQ, BestQ = self.elements._AvQ_, \
                             self.elements._WorstQ_, \
                             self.elements._BestQ_

        AvQ = cOmm.gather(AvQ, root=mAster_rank)
        WorstQ = cOmm.gather(WorstQ, root=mAster_rank)
        BestQ = cOmm.gather(BestQ, root=mAster_rank)

        if rAnk == mAster_rank:
            AvQ = [i for i in AvQ if i]
            WorstQ = [i for i in WorstQ if i]
            BestQ = [i for i in BestQ if i]
            AvQ = sum(AvQ) / len(AvQ)
            WorstQ = min(WorstQ)
            BestQ = max(BestQ)
        else:
            pass
        AvQ = cOmm.bcast(AvQ, root=mAster_rank)
        WorstQ = cOmm.bcast(WorstQ, root=mAster_rank)
        BestQ = cOmm.bcast(BestQ, root=mAster_rank)

        return {'average quality': AvQ,
                'worst quality': WorstQ,
                'best quality': BestQ}


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
        self.___PRIVATE_generating_trace_elements___()
        self.RESET_cache()
        self._freeze_self_()


    def RESET_cache(self):
        self._type_amount_dict_ = None
        self._Q_ = None
        self._AvQ_ = None
        self._WorstQ_ = None
        self._BestQ_ = None
        self.___multi_elements_metric___ = None


    def ___DO_find_type_and_amount_numbered_before___(self):
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
    def DO(self):
        if self._DO_ is None:
            self._DO_ = _3dCSCG_Trace_Elements_DO(self)
        return self._DO_

    @property
    def SELFCHECK(self):
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
                self.___PRIVATE_get_quality_of_trace_elements___
        return self._Q_

    @property
    def coordinate_transformation(self):
        if self._ct_ is None:
            self._ct_ = _3dCSCG_Trace_Elements_CoordinateTransformation(self)
        return self._ct_

    @property
    def ___PRIVATE_get_quality_of_trace_elements___(self):
        """

        :return: A dict whose keys are the numbers of local trace
            elements and values are the qualities of the trace elements.
        """
        r = np.array([-0.5, -0.5, 0.5, 0.5, 0, 0.75, 0.75, -0.75, -0.75])
        s = np.array([-0.5, 0.5, -0.5, 0.5, 0, -0.75, 0.75, -0.75, 0.75])

        xi, et, sg = np.array([0, ]), np.array([0, ]), np.array([0, ])
        rc, sc = np.array([0, ]), np.array([0, ])

        LEN = len(r)
        assert len(s) == LEN, f"r, s length dis-match."
        Q = dict()
        Qs = list()
        for i in self:
            te = self[i]
            e = te.CHARACTERISTIC_element
            side = te.CHARACTERISTIC_side

            uv = te.coordinate_transformation.___PRIVATE_outward_unit_normal_vector___(r, s, from_element=e, side=side)
            vx, vy, vz = uv
            V = [[vx[_], vy[_], vz[_]] for _ in range(LEN)]
            comb = combinations(V, 2)
            angle = list()
            for v1v2 in comb:
                angle.append(angle_between_two_vectors(*v1v2))
            angle = max(angle)
            quality_0 =  float(1 - angle / np.pi)

            me = self._mesh_.elements[e]
            me_ct = me.coordinate_transformation.mapping(xi, et, sg)
            te_ct = te.coordinate_transformation.mapping(rc, sc, from_element=e, side=side)
            te_uv = te.coordinate_transformation.___PRIVATE_outward_unit_normal_vector___(rc, sc, from_element=e, side=side)
            v_inner = (me_ct[0]-te_ct[0], me_ct[1]-te_ct[1], me_ct[2]-te_ct[2])
            v_outer = te_uv
            angle = angle_between_two_vectors(v_inner, v_outer)
            quality_1 = float(angle / np.pi)


            _ = min([quality_0, quality_1])
            Q[i] = _
            Qs.append(_)

        if len(Qs) > 0:
            AvQ = sum(Qs) / len(Qs)
            WorstQ = min(Qs)
            BestQ = max(Qs)
        else:
            AvQ = None
            WorstQ = None
            BestQ = None

        return Q, AvQ, WorstQ, BestQ

    @property
    def _multi_elements_metric_(self):
        """"""
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
        (int) The total number of trace elements.

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
            between [-1,1], we will pick up two from them according the the side and do the mesh grid.

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

    def ___DO_illustrate_trace_element___(self, i, density_factor=2):
        """We illustrate the trace element #i with matplotlib on the
        mesh elements which share this trace element.

        To call this method, call it from all cores, otherwise, it does
        not work.

        :param int i: The trace element #i will be illustrated.
        :param int density_factor: be in [1,10], not be too large.
        :return:
        """
        # wo first find which core is going to do the plot _____________
        if i in self:
            in_or_out = rAnk
            I_am_in = True
        else:
            in_or_out = -1
            I_am_in = False

        who_are_in = cOmm.gather(in_or_out, root=mAster_rank)
        if rAnk != mAster_rank:
            who_will_do_it = None
            in_how_many_cores = None
            the_other_core = None
        else:
            who_will_do_it = max(who_are_in)
            in_how_many_cores = 0
            the_other_core = None
            for _ in who_are_in:
                if _ != -1:
                    if _ !=who_will_do_it:
                        assert the_other_core is None
                        the_other_core = _
                    in_how_many_cores += 1

        who_will_do_it = cOmm.bcast(who_will_do_it, root=mAster_rank)
        in_how_many_cores = cOmm.bcast(in_how_many_cores, root=mAster_rank)
        the_other_core = cOmm.bcast(the_other_core, root=mAster_rank)

        sent_to = recv_from = None

        if in_how_many_cores == 1:
            if rAnk == who_will_do_it:
                assert self[i].IS_shared_by_cores is False, "trivial check"
            else:
                assert I_am_in is False, "trivial check"
        elif in_how_many_cores == 2:
            if I_am_in and rAnk != who_will_do_it:
                assert rAnk == the_other_core, "trivial check"
                # print(rAnk, ': I am going to provide data to rAnk ', who_will_do_it)
                sent_to = who_will_do_it
            elif I_am_in and rAnk == who_will_do_it:
                recv_from = the_other_core
            else:
                I_am_in = False
        else:
            raise Exception(f"can only have two cores sharing one trace element!")

        DATA2 = None
        r = s = cs = cr = uv_r = uv_s = None
        # prepare coordinate data first _____________________________________________
        if rAnk == who_will_do_it or I_am_in:
            density = 5 + 4 * density_factor
            i0 = 1 + density_factor
            i1 = 2 * density_factor + 2
            i2 = 3 * density_factor + 3
            _ = np.linspace(-1, 1, density)
            r, s = np.meshgrid(_, _, indexing='ij')
            anchors = ( # the points we will plot the outward unit norm vector
                [i0, i0],
                [i0, i2],
                [i2, i0],
                [i2, i2],
                [i1, i1],
            )
            uv_r = np.array([r[indices[0],indices[1]] for indices in anchors])
            uv_s = np.array([s[indices[0],indices[1]] for indices in anchors])
            cr, cs = np.array([r[i1, i1],]), np.array([s[i1, i1],])
            np.testing.assert_almost_equal(cr, 0)
            np.testing.assert_almost_equal(cs, 0)

        # Now lets prepare the data for the second subplot______________
        if in_how_many_cores == 1:
            if rAnk == who_will_do_it:
                if self[i].IS_on_mesh_boundary:
                    DATA2 = self[i].NON_CHARACTERISTIC_position # a string
                else:
                    ncp = self[i].NON_CHARACTERISTIC_position
                    other_element = int(ncp[:-1])
                    element = other_element
                    side = ncp[-1]
                    tes = self.map[element]
                    DATA2 = dict()
                    DATA2['element'] = element
                    DATA2['side'] = side
                    DATA2['tes'] = tes
                    for _, tei in enumerate(tes):
                        _side_ = 'NSWEBF'[_]
                        te = self._elements_[tei]
                        X, Y, Z = te.coordinate_transformation.mapping(r, s, from_element=element, side=_side_)

                        if tei == i and _side_ == side:
                            x, y, z = te.coordinate_transformation.mapping(uv_r, uv_s, from_element=element, side=_side_)
                            uvx, uvy, uvz = te.coordinate_transformation.___PRIVATE_outward_unit_normal_vector___(uv_r, uv_s, from_element=element, side=_side_)
                            DATA2[str(tei)+_side_] = [(X, Y, Z), (x, y, z), (uvx, uvy, uvz)]
                        else:
                            x, y, z = te.coordinate_transformation.mapping(cr, cs, from_element=element, side=_side_)
                            DATA2[str(tei)+_side_] = [(X, Y, Z), (x, y, z), _side_]

            else:
                pass
        else: # in_how_many_cores == 2
            if rAnk == the_other_core:
                element = self[i].CHARACTERISTIC_element
                side = self[i].CHARACTERISTIC_side
                tes = self.map[element]
                DATA2 = dict()
                DATA2['element'] = element
                DATA2['side'] = side
                DATA2['tes'] = tes
                # noinspection PyTypeChecker
                assert tes.count(i) == 1, f"must be the case."
                for _, tei in enumerate(tes):
                    _side_ = 'NSWEBF'[_]

                    te = self._elements_[tei]
                    X, Y, Z = te.coordinate_transformation.mapping(r, s, from_element=element, side=_side_)
                    if tei == i:
                        x, y, z = te.coordinate_transformation.mapping(uv_r, uv_s, from_element=element, side=_side_)
                        uvx, uvy, uvz = te.coordinate_transformation.___PRIVATE_outward_unit_normal_vector___(uv_r, uv_s, from_element=element, side=_side_)
                        DATA2[str(tei)+_side_] = [(X, Y, Z), (x, y, z), (uvx, uvy, uvz)]
                    else:
                        x, y, z = te.coordinate_transformation.mapping(cr, cs, from_element=element, side=_side_)
                        DATA2[str(tei)+_side_] = [(X, Y, Z), (x, y, z),  _side_   ]

                DATA2['rank'] = rAnk
                cOmm.send(DATA2, dest=sent_to, tag=who_will_do_it)
            elif rAnk == who_will_do_it:
                DATA2 = cOmm.recv(source=recv_from, tag=who_will_do_it)
            else:
                pass

        # lets do the plots ____________________________________________
        if rAnk == who_will_do_it:
            element = self[i].CHARACTERISTIC_element
            side = self[i].CHARACTERISTIC_side
            tes = self.map[element]
            assert i in tes, "trivial check"
            # noinspection PyTypeChecker
            assert 'NSWEBF'[tes.index(i)] == side

            fig = plt.figure(figsize=(14,6))
            TITLE = f"**trace element {i}**"
            if self[i].IS_on_periodic_boundary:
                TITLE += ' (periodic)'
            fig.suptitle(TITLE)

            ax = fig.add_subplot(121, projection='3d')
            for _, tei in enumerate(tes):
                te = self._elements_[tei]
                _side_ = 'NSWEBF'[_]
                X, Y, Z = te.coordinate_transformation.mapping(r, s, from_element=element, side=_side_)
                if tei == i and _side_ == side:
                    ax.plot_surface(X, Y, Z) # plot the trace element
                    x, y, z = te.coordinate_transformation.mapping(uv_r, uv_s, from_element=element, side=_side_)
                    uvx, uvy, uvz = te.coordinate_transformation.___PRIVATE_outward_unit_normal_vector___(uv_r, uv_s, from_element=element, side=_side_)

                    x_range = np.max(X) - np.min(X)
                    y_range = np.max(Y) - np.min(Y)
                    z_range = np.max(Z) - np.min(Z)
                    mean_range = (x_range + y_range + z_range) / 6

                    ax.quiver(x, y, z, uvx*mean_range, uvy*mean_range, uvz*mean_range, color='r', linewidth=0.5)
                    ax.text(x[-1] + 0.5*uvx[-1]*mean_range, y[-1] + 0.5*uvy[-1]*mean_range, z[-1] + 0.5*uvz[-1]*mean_range,
                            side, color='green', ha='center', va='center', ma='center')
                else:
                    ax.plot_surface(X, Y, Z, color=(0.7,0.7,0.7,0.5))
                    x, y, z = te.coordinate_transformation.mapping(cr, cs, from_element=element, side=_side_)
                    ax.text(x[0], y[0], z[0], 'NSWEBF'[_],
                            color='k', ha='center', va='center', ma='center')

            ax.set_xlabel(r'$x$')
            ax.set_ylabel(r'$y$')
            ax.set_zlabel(r'$z$')

            if in_how_many_cores == 1 and isinstance(DATA2, str):
                plt.title(f"in rank {rAnk}, on [{side}] of element {element}.")
            elif in_how_many_cores == 1 and isinstance(DATA2, dict):
                plt.title(f"in rank {rAnk}, on [{side}] of characteristic element {element}.")
            elif in_how_many_cores == 2:
                plt.title(f"in rank {rAnk}, on [{side}] of characteristic element {element}.")
            else:
                raise Exception()

            ax = fig.add_subplot(122, projection='3d')
            if isinstance(DATA2, dict):
                element = DATA2['element']
                side_0 = side
                side = DATA2['side']
                tes = DATA2['tes']
                assert i in tes, "trivial check"
                # noinspection PyTypeChecker
                # assert 'NSWEBF'[tes.index(i)] == side
                assert side+side_0 in ('NS', 'SN', 'WE', 'EW', 'BF', 'FB')
                for _, tei in enumerate(tes):
                    _side_ = 'NSWEBF'[_]
                    if tei == i and _side_ == side:
                        XYZ, xyz, uv_xyz = DATA2[str(tei)+_side_]
                        ax.plot_surface(*XYZ)  # plot the trace element
                        uvx, uvy, uvz = uv_xyz
                        x, y, z = xyz
                        X, Y, Z = XYZ
                        x_range = np.max(X) - np.min(X)
                        y_range = np.max(Y) - np.min(Y)
                        z_range = np.max(Z) - np.min(Z)
                        mean_range = (x_range + y_range + z_range) / 6
                        ax.quiver(*xyz, uvx*mean_range, uvy*mean_range, uvz*mean_range, color='r', linewidth=0.5)
                        ax.text(x[-1] + 0.5*uvx[-1]*mean_range, y[-1] + 0.5*uvy[-1]*mean_range, z[-1] + 0.5*uvz[-1]*mean_range,
                                side, color='green', ha='center', va='center', ma='center')
                    else:
                        XYZ, xyz, _s_ = DATA2[str(tei)+_side_]
                        ax.plot_surface(*XYZ, color=(0.7, 0.7, 0.7, 0.5))
                        x, y, z = xyz
                        ax.text(x[0], y[0], z[0], _s_,
                                color='k', ha='center', va='center', ma='center')

            ax.set_xlabel(r'$x$')
            ax.set_ylabel(r'$y$')
            ax.set_zlabel(r'$z$')
            if in_how_many_cores == 1 and isinstance(DATA2, str):
                plt.title(f"on mesh boundary: <{DATA2}>")
            elif in_how_many_cores == 1 and isinstance(DATA2, dict):
                plt.title(f"in rank {rAnk}, on [{side}] of element {element}.")
            elif in_how_many_cores == 2:
                plt.title(f"in rank {DATA2['rank']}, on [{side}] of characteristic element {element}.")
            else:
                raise Exception()
            plt.show()

    def ___SELFCHECK_outward_unit_normal_vector___(self):
        """This is a self check program to check that we obtain the
        correct outward unit normal vector(s) for a trace element.

        Notice this self check program does not give error. It will only
        give warning; it can only find something may be wrong.

        We basically compute the angle between the outward norm vector
        and the vector pointing the mesh element center. If the angle
        is pi (180 degree), of course, the outward normal vector is
        correct. If it is 0, then we have a wrong outward normal vector.
        """
        # the center of the mesh element.
        xi, et, sg = np.array([0,]), np.array([0,]), np.array([0,])
        # the center of the trace element.
        rc, sc = np.array([0,]), np.array([0,])
        for i in self:
            te = self[i]
            e1 = te.CHARACTERISTIC_element
            s1 = te.CHARACTERISTIC_side
            E = [e1,]
            S = [s1,]
            if not te.IS_on_mesh_boundary:
                p2 = te.NON_CHARACTERISTIC_position
                e2 = int(p2[:-1])
                s2 = p2[-1]
                if e2 in self._mesh_.elements:
                    E.append(e2)
                    S.append(s2)
            for e, s in zip(E, S):
                me = self._mesh_.elements[e]
                me_ct = me.coordinate_transformation.mapping(xi, et, sg)
                te_ct = te.coordinate_transformation.mapping(rc, sc, from_element=e, side=s)
                te_uv = te.coordinate_transformation.___PRIVATE_outward_unit_normal_vector___(rc, sc, from_element=e, side=s)

                v_inner = (me_ct[0]-te_ct[0], me_ct[1]-te_ct[1], me_ct[2]-te_ct[2])
                v_outer = te_uv
                angle = angle_between_two_vectors(v_inner, v_outer)

                if angle < np.pi/4:
                    warnings.warn(
                        f"___PRIVATE_outward_unit_normal_vector___ of trace element "
                        f"#{i} may be wrong",
                        TraceElementWarning)
                elif angle < np.pi/2:
                    warnings.warn(
                        f"the mesh element #{e} may be very distorted.",
                        TraceElementWarning)
                else:
                    pass


class _3dCSCG_Trace_Elements_DO(FrozenOnly):
    """We find some specific groups of elements."""

    def __init__(self, trace_elements):
        self._elements_ = trace_elements
        self._freeze_self_()

    def illustrate_trace_element(self, *args, **kwargs):
        return self._elements_.___DO_illustrate_trace_element___(*args, **kwargs)


class _3dCSCG_Trace_Elements_SELFCHECK(FrozenOnly):
    """We find some specific groups of elements."""

    def __init__(self, trace_elements):
        self._elements_ = trace_elements
        self._freeze_self_()

    def outward_unit_normal_vector(self, *args, **kwargs):
        return self._elements_.___SELFCHECK_outward_unit_normal_vector___(*args, **kwargs)


class _3dCSCG_Trace_Elements_Group(FrozenOnly):
    """We find some specific groups of elements."""
    def __init__(self, trace_elements):
        self._elements_ = trace_elements
        self._freeze_self_()

    def elements_on_same_plane_as(self, i):
        """We find all trace elements are topologically on the same
        plane as trace element No. ``i``.

        For example, in the crazy mesh of 2*2*2 mesh elements. Along
        the x-axis, there are 12 trace elements in total  perpendicular
        to x-axis. They are distributed on 3 levels; 4 trace elements on
        each level. If on the middle level, the 4 trace elements are
        numbered 1, 4, 7, 10. Now we do
        ``mesh.trace.elements.group.elements_on_same_plane_as(i)``,
        where i is in {1,4,7,10}, we get {1,4,7,10}. So, we get all
        trace elements that are on the same plane topologically.

        :param int i:
        """
        raise NotImplementedError()


class _3dCSCG_Trace_Elements_CoordinateTransformation(FrozenOnly):
    def __init__(self, trace_elements):
        self._elements_ = trace_elements
        self._freeze_self_()

    def mapping(self, *args, **kwargs):
        """"""
        raise Exception("Can not access mapping from trace elements, must do it from particular trace element.")

    def Jacobian_matrix(self, *evaluation_points_3):
        """
        The local Jacobian matrix.

        :param evaluation_points_3: A tuple or list of shape (ndim, ...).
        """
        return TraceElementsCTValuesCache(self._elements_, 'Jacobian_matrix', evaluation_points_3)

    def inverse_Jacobian_matrix(self, *evaluation_points_3):
        """
        The local inverse_Jacobian matrix.

        :param evaluation_points_3 : A tuple or list of shape (ndim, ...).
        """
        return TraceElementsCTValuesCache(self._elements_, 'inverse_Jacobian_matrix', evaluation_points_3)

    def metric_matrix(self, *evaluation_points_3):
        """Compute the metric matrix G."""
        return TraceElementsCTValuesCache(self._elements_, 'metric_matrix', evaluation_points_3)

    def metric(self, *evaluation_points_3):
        """return metric g."""
        return TraceElementsCTValuesCache(self._elements_, 'metric', evaluation_points_3)

    def unit_normal_vector(self, *evaluation_points_3):
        """return metric g."""
        return TraceElementsCTValuesCache(self._elements_, 'unit_normal_vector', evaluation_points_3)

class TraceElementsCTValuesCache(FrozenOnly):
    def __init__(self, trace_elements, CTT, evaluation_points_3):
        self._elements_ = trace_elements
        assert isinstance(CTT, str) and CTT != 'mapping', f"CTT={CTT} wrong."
        self._CTT_ = CTT
        self._0ep_ = [evaluation_points_3[1], evaluation_points_3[2]]
        self._e1p_ = [evaluation_points_3[0], evaluation_points_3[2]] # do not use 2, 0
        self._ep2_ = [evaluation_points_3[0], evaluation_points_3[1]]
        self._multi_elements_metric_ = self._elements_._multi_elements_metric_
        self._cache_ = dict()
        self._freeze_self_()

    def __getitem__(self, i):
        element = self._elements_[i]
        type_wrt_metric = element.type_wrt_metric.mark

        if isinstance(type_wrt_metric, int): # it is unique, we use the id (int) as the mark. Otherwise, it must be a string.
            side = element.CHARACTERISTIC_side
            if side in 'NS':
                _xi_eta_sigma_ = self._0ep_
            elif side in 'WE':
                _xi_eta_sigma_ = self._e1p_
            elif side in 'BF':
                _xi_eta_sigma_ = self._ep2_
            else:
                raise Exception()
            return getattr(element.coordinate_transformation, self._CTT_)(
                *_xi_eta_sigma_)

        elif type_wrt_metric in self._cache_:
            return self._cache_[type_wrt_metric]

        else:

            side = element.CHARACTERISTIC_side
            if side in 'NS':
                _xi_eta_sigma_ = self._0ep_
            elif side in 'WE':
                _xi_eta_sigma_ = self._e1p_
            elif side in 'BF':
                _xi_eta_sigma_ = self._ep2_
            else:
                raise Exception()
            result = getattr(element.coordinate_transformation, self._CTT_)(
                *_xi_eta_sigma_)

            if type_wrt_metric in self._multi_elements_metric_ and \
                type_wrt_metric not in self._cache_ and \
                self._multi_elements_metric_[type_wrt_metric] >= caChe_factor:
                # here we have very strict cache rule.
                self._cache_[type_wrt_metric] = result

            return result

    def __len__(self):
        return len(self._elements_)

    def __contains__(self, item):
        return item in self._elements_

    def __iter__(self):
        for i in self._elements_:
            yield i












class _3dCSCG_Trace_Element(FrozenOnly):
    """

    :param trace_elements:
    :param i:
    :param position_1:
    :param position_2:
    :param cp: characteristic position
    :param ondb:
    :param onpb:
    """
    def __init__(self, trace_elements, i, position_1, position_2, cp, ondb=False, onpb=False):
        self._elements_ = trace_elements
        self._mesh_ = trace_elements._mesh_
        self._i_ = i
        self._p1_ = position_1
        self._p2_ = position_2
        self._cp_ = cp
        if position_1 == cp:
            self._ncp_ = position_2
        elif position_2 == cp:
            self._ncp_ = position_1
        else:
            raise Exception()
        self._ondb_ = ondb
        self._onpb_ = onpb
        assert self.CHARACTERISTIC_element in self._elements_._mesh_.elements, \
            "CHARACTERISTIC_element must be in the same core."
        if self.IS_on_mesh_boundary:
            assert self.NON_CHARACTERISTIC_position[0] not in '1234567890'
        self._ct_ = None
        self._type_wrt_metric_ = None
        self._freeze_self_()
        # # do a check for periodic trace element ________________________
        # if self.IS_on_periodic_boundary:
        #     assert not self.IS_on_mesh_boundary # must be the case
        #     e1 = int(position_1[:-1])
        #     e2 = int(position_2[:-1])
        #     if e1 == e2:
        #         warnings.warn(f"periodic trace element #{i} is on two "
        #                       f"faces of same mesh element #{e1}, "
        #                       f"this may cause unknown error. Please "
        #                       f"consider use more mesh elements to "
        #                       f"remove this situation.",
        #                       PeriodicTraceElementWarning)

    @property
    def positions(self):
        return self._p1_, self._p2_

    @property
    def coordinate_transformation(self):
        if self._ct_ is None:
            self._ct_ = _3dCSCG_Trace_Element_CoordinateTransformation(self)
        return self._ct_

    @property
    def normal_direction(self):
        """"""
        if self._p1_[-1] in 'NS':
            return 'NS'
        elif self._p1_[-1] in 'WE':
            return 'WE'
        elif self._p1_[-1] in 'BF':
            return 'BF'
        else:
            raise Exception()

    @property
    def NON_CHARACTERISTIC_position(self):
        """The other position."""
        return self._ncp_

    @property
    def CHARACTERISTIC_position(self):
        """The position we mainly locate this trace element."""
        return self._cp_
    @property
    def CHARACTERISTIC_element(self):
        """We mainly consider this trace element is a side of this mesh
        element."""
        return int(self._cp_[:-1])
    @property
    def CHARACTERISTIC_side(self):
        """We main consider this trace element is such a side of the
        CHARACTERISTIC_element."""
        return self._cp_[-1]
    @property
    def CHARACTERISTIC_region(self):
        """We mainly consider this trace element is in this region."""
        region = self._mesh_.DO.FIND_region_name_of_element(
            self.CHARACTERISTIC_element)
        return region

    @property
    def spacing(self):
        element_spacing = self._mesh_.elements[self.CHARACTERISTIC_element].spacing
        side = self.CHARACTERISTIC_side
        if side == 'N':
            trace_spacing = (element_spacing[0][0], element_spacing[1], element_spacing[2])
        elif side == 'S':
            trace_spacing = (element_spacing[0][1], element_spacing[1], element_spacing[2])
        elif side == 'W':
            trace_spacing = (element_spacing[0], element_spacing[1][0], element_spacing[2])
        elif side == 'E':
            trace_spacing = (element_spacing[0], element_spacing[1][1], element_spacing[2])
        elif side == 'B':
            trace_spacing = (element_spacing[0], element_spacing[1], element_spacing[2][0])
        elif side == 'F':
            trace_spacing = (element_spacing[0], element_spacing[1], element_spacing[2][1])
        else:
            raise Exception()
        return trace_spacing

    @property
    def i(self):
        """This is the ith trace element."""
        return self._i_

    @property
    def IS_on_mesh_boundary(self):
        """As this property name says."""
        return self._ondb_

    @property
    def on_mesh_boundary(self):
        """Return the mesh boundary name this trace element is on. If it is not on one, return None."""
        if self.IS_on_mesh_boundary:
            return self.NON_CHARACTERISTIC_position
        else:
            return None

    @property
    def IS_on_periodic_boundary(self):
        """As this property name says."""
        return self._onpb_

    @property
    def IS_shared_by_cores(self):
        """True or False, As this property name says."""
        if self.IS_on_mesh_boundary:
            return False
        else:
            if int(self._p1_[:-1]) in self._elements_._mesh_.elements and \
                int(self._p2_[:-1]) in self._elements_._mesh_.elements:
                return False
            else:
                return True

    @property
    def shared_with_core(self):
        """None or a int. """
        if self.IS_shared_by_cores:
            if int(self._p1_[:-1]) in self._elements_._mesh_.elements:
                CORE = self._elements_._mesh_.DO.FIND_slave_of_element(int(self._p2_[:-1]))
            elif int(self._p2_[:-1]) in self._elements_._mesh_.elements:
                CORE = self._elements_._mesh_.DO.FIND_slave_of_element(int(self._p1_[:-1]))
            else:
                raise Exception()
            assert CORE != rAnk
            return CORE
        else:
            return None

    @property
    def type_wrt_metric(self):
        """Return the trace-element-metric-type object reflecting the element type."""
        if self._type_wrt_metric_ is None:

            self._type_wrt_metric_ = \
                self._mesh_.domain.regions[
                    self.CHARACTERISTIC_region].type_wrt_metric.___CLASSIFY_TRACE_ELEMENT_of_spacing___(
                    self.spacing)

        return self._type_wrt_metric_






class _3dCSCG_Trace_Element_CoordinateTransformation(FrozenOnly):
    def __init__(self, te):
        self._te_ = te
        self._freeze_self_()

    def mapping(self, *evaluation_points, from_element=None, side=None, parse_3_1d_eps=False):
        """
        The local mapping.

        :param evaluation_points : A tuple or list of shape (ndim-1, ...).
        :param from_element: (default: ``None``) We try to compute the
            mapping from a given mesh element. When it is None,
            we compute it from its CHARACTERISTIC_element. Notice that
            when a trace element in on periodic domain, we have to
            provide a ``from_element``. Because it this is the case,
            we will get different mapping from different element.
        :param side: when one mesh element is periodic to itself, we may
            need to provide side as well.
        :param parse_3_1d_eps: If `parse_ep` is True, then we have *ep is xi, eta, sigma, and they are all 1d,
            between [-1,1], we will pick up two from them according the the side and do the mesh grid.
        """
        if self._te_.IS_on_periodic_boundary:
            assert from_element is not None, \
                "to compute the mapping of a trace element on periodic " \
                "boundary, we have to provide from which element you " \
                "want to compute it since it clearly will gives " \
                "different results."
            if from_element == 'any':
                # the different results do not matter; for example, when
                # we want to get value from a periodic function, the
                # location for evaluating the function also does not
                # matter.
                from_element = self._te_.CHARACTERISTIC_element
            else:
                pass

        if from_element is None:
            i = self._te_.CHARACTERISTIC_element
        elif from_element == 'any':
            i = self._te_.CHARACTERISTIC_element
        else:
            i = from_element

        assert self._te_.i in self._te_._elements_.map[i], \
            f"trace element {self._te_.i} is not on mesh element {i}."

        if self._te_._elements_.map[i].count(self._te_.i) == 1:
            side_index = self._te_._elements_.map[i].index(self._te_.i)
            element_side = 'NSWEBF'[side_index]
            if from_element is None:
                assert element_side == self._te_.CHARACTERISTIC_side

            if side is not None:
                assert element_side == side, f"cannot compute on provided side {side}"

        elif self._te_._elements_.map[i].count(self._te_.i) == 2:
            assert side is not None, f"trace element #{self._te_.i} " \
                                     f"is on two sides of element #{i} " \
                                     f"(periodic), provide side as well."
            element_side = side
        else:
            raise Exception()

        assert self._te_.i == self._te_._elements_.map[i]['NSWEBF'.index(element_side)], \
            f"trace element #{self._te_.i} is not at {element_side} of mesh element #{i}."

        ep = self._te_._elements_.___generate_full_ep___(evaluation_points, element_side, parse_3_1d_eps=parse_3_1d_eps)
        x, y, z = self._te_._mesh_.elements[i].coordinate_transformation.mapping(*ep)
        return x, y, z

    def Jacobian_matrix(self, *evaluation_points, parse_3_1d_eps=False):
        """
        The local Jacobian matrix.

        :param evaluation_points: A tuple or list of shape (ndim-1, ...).
        :param parse_3_1d_eps:
        """
        i = self._te_.CHARACTERISTIC_element
        element_side = self._te_.CHARACTERISTIC_side
        ep = self._te_._elements_.___generate_full_ep___(evaluation_points, element_side, parse_3_1d_eps=parse_3_1d_eps)
        J = self._te_._mesh_.elements[i].coordinate_transformation.Jacobian_matrix(*ep)
        if element_side in 'NS':
            return ((J[0][1], J[0][2]),
                    (J[1][1], J[1][2]),
                    (J[2][1], J[2][2]))
        elif element_side in 'WE': # this is very important, do not used [0], [2] for the indices. But when we feed x, y, z, always use x, y, z.
            return ((J[0][2], J[0][0]),
                    (J[1][2], J[1][0]),
                    (J[2][2], J[2][0]))
        else:
            return ((J[0][0], J[0][1]),
                    (J[1][0], J[1][1]),
                    (J[2][0], J[2][1]))

    def inverse_Jacobian_matrix(self, *evaluation_points):
        """
        The local inverse_Jacobian matrix.

        :param evaluation_points : A tuple or list of shape (ndim-1, ...).
        """
        i = self._te_.CHARACTERISTIC_element
        element_side = self._te_.CHARACTERISTIC_side
        ep = self._te_._elements_.___generate_full_ep___(evaluation_points, element_side)
        iJ = self._te_._mesh_.elements[i].coordinate_transformation.inverse_Jacobian_matrix(*ep)
        if element_side in 'NS':
            return ((iJ[1][0], iJ[1][1], iJ[1][2]),
                    (iJ[2][0], iJ[2][1], iJ[2][2]))
        elif element_side in 'WE': # this is very important, do not used [0], [2] for the indices.. But when we feed x, y, z, always use x, y, z.
            return ((iJ[2][0], iJ[2][1], iJ[2][2]),
                    (iJ[0][0], iJ[0][1], iJ[0][2]))
        else:
            return ((iJ[0][0], iJ[0][1], iJ[0][2]),
                    (iJ[1][0], iJ[1][1], iJ[1][2]))

    def metric_matrix(self, *evaluation_points):
        """Compute the metric matrix G whose entries are g_{i,j}."""
        J = self.Jacobian_matrix(*evaluation_points)
        Gk = [[None for _ in range(2)] for __ in range(2)]
        for i in range(2):
            for j in range(i, 2):
                Gk[i][j] = J[0][i] * J[0][j]
                for l in range(1, 3):
                    Gk[i][j] += J[l][i] * J[l][j]
                if i != j:
                    Gk[j][i] = Gk[i][j]
        return Gk


    def inverse_metric_matrix(self, *evaluation_points):
        """Compute the inverse metric matrix G whose entries are g^{i,j}."""
        iJ = self.inverse_Jacobian_matrix(*evaluation_points)
        Gk = [[None for _ in range(2)] for __ in range(2)]
        for i in range(2):
            for j in range(i, 2):
                Gk[i][j] = iJ[i][0] * iJ[j][0] + iJ[i][1] * iJ[j][1] + iJ[i][2] * iJ[j][2]
                if i != j:
                    Gk[j][i] = Gk[i][j]

        return Gk



    def metric(self, *evaluation_points):
        """return metric g."""
        G = self.metric_matrix(*evaluation_points)
        # noinspection PyUnresolvedReferences
        return G[0][0]*G[1][1] - G[0][1]*G[1][0]

    def ___PRIVATE_outward_unit_normal_vector___(self, *evaluation_points, from_element=None, side=None):
        """For a single trace element (a surface in 3D space), it is
        hard to say which direction is the positive direction (the
        outward direction). So we put onto the surface of the mesh
        element. So, this outward unit norm vector of this trace element
        is the outward unit norm vector of the characteristic mesh
        element on the face the trace element is located.

        In fact, for a trace element, it seems not good to consider it
        has a normal vector. We should consider so for a mesh element.
        For example, when we have a 2-form vector(u), the trace of it,
        vector(u) dot vector(n), gives a scalar field anyway, the
        normal vector vector(n) is used at the mesh element level, as
        well as the trace matrix N.

        vec(a) = (a1, a2, a3)
        vec(b) = (b1, b2, b3)

        vec(a) x vec(b) = (a2 b3 - a3 b2, a3 b1 - a1 b3, a1 b2-a2 b1)

        :param evaluation_points: A tuple or list of shape (ndim-1, ...).
        :param from_element: We will compute the normal vector by
            considering it is on this mesh element. Please give this
            parameter in case of random error!
        :param side: when one mesh element is periodic to itself, we may
            need to provide side as well.
        """
        J = self.Jacobian_matrix(*evaluation_points)
        a = (J[0][0], J[1][0], J[2][0])
        b = (J[0][1], J[1][1], J[2][1])
        acb0 = a[1] * b[2] - a[2] * b[1]
        acb1 = a[2] * b[0] - a[0] * b[2]
        acb2 = a[0] * b[1] - a[1] * b[0]
        norm = np.sqrt(acb0**2 + acb1**2 + acb2**2)

        if from_element is None:
            i = self._te_.CHARACTERISTIC_element
        else:
            i = from_element
        assert self._te_.i in self._te_._elements_.map[i], \
            f"trace element {self._te_.i} is not on mesh element {i}."

        if self._te_._elements_.map[i].count(self._te_.i) == 1:
            side_index = self._te_._elements_.map[i].index(self._te_.i)
            element_side = 'NSWEBF'[side_index]
            if from_element is None:
                assert element_side == self._te_.CHARACTERISTIC_side

            if side is not None:
                assert element_side == side, f"cannot compute on provided side {side}"

        elif self._te_._elements_.map[i].count(self._te_.i) == 2:
            assert side is not None, f"trace element #{self._te_.i} " \
                                     f"is on two sides of element #{i} " \
                                     f"(periodic), provide side as well."
            element_side = side
        else:
            raise Exception()

        if side is not None: element_side = side #

        uv = np.array([acb0 / norm, acb1 / norm, acb2 / norm])
        if element_side in 'NWB':
            uv = - uv
        else:
            pass

        return uv

    def unit_normal_vector(self, *evaluation_points, parse_3_1d_eps=False):
        """The unit_normal_vector (vector n) according to the right-hand-rule (without considering it is
        which side of the mesh element; it can point inner direction of the mesh element).

        vec(a) = (a1, a2, a3)
        vec(b) = (b1, b2, b3)

        vec(a) x vec(b) = (a2 b3 - a3 b2, a3 b1 - a1 b3, a1 b2-a2 b1)

        :param evaluation_points: A tuple or list of shape (ndim-1, ...).
        :param parse_3_1d_eps:
        """
        J = self.Jacobian_matrix(*evaluation_points, parse_3_1d_eps=parse_3_1d_eps)
        a = (J[0][0], J[1][0], J[2][0])
        b = (J[0][1], J[1][1], J[2][1])
        acb0 = a[1] * b[2] - a[2] * b[1]
        acb1 = a[2] * b[0] - a[0] * b[2]
        acb2 = a[0] * b[1] - a[1] * b[0]
        norm = np.sqrt(acb0**2 + acb1**2 + acb2**2)

        return np.array([acb0 / norm, acb1 / norm, acb2 / norm])








if __name__ == '__main__':
    # mpiexec -n 12 python _3dCSCG\mesh\trace.py
    from _3dCSCG.main import MeshGenerator
    elements = [3, 4, 2]
    mesh = MeshGenerator('crazy_periodic', c=0.3, bounds=([0,1], [0,1], [0,1]))(elements)
    mesh.trace.elements.SELFCHECK.outward_unit_normal_vector()
    Q = mesh.trace.elements.quality
    print(mesh.quality)
    print(mesh.trace.quality)

    mesh.trace.elements.DO.illustrate_trace_element(1)

    # te0 = mesh.trace.elements[0]

    # print(te0.IS_on_periodic_boundary)

    for i in range(mesh.trace.elements.GLOBAL_num):
        mesh.trace.elements.DO.illustrate_trace_element(i)
        # if i in mesh.trace.elements:
        #     te = mesh.trace.elements[i]
        #
        #     print(rAnk, te.type_wrt_metric.mark)