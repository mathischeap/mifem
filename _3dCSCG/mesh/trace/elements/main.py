



import sys
if './' not in sys.path: sys.path.append('./')

from root.config import *
from SCREWS.frozen import FrozenOnly



from SCREWS.customized_warnings import TraceElementWarning
from SCREWS.functions._3d import angle_between_two_vectors
from itertools import combinations

from _3dCSCG.mesh.trace.elements.DO import _3dCSCG_Trace_Elements_DO
from _3dCSCG.mesh.trace.elements.group import _3dCSCG_Trace_Elements_Group
from _3dCSCG.mesh.trace.elements.coordinate_transformation import _3dCSCG_Trace_Elements_CoordinateTransformation
from _3dCSCG.mesh.trace.elements.element.main import _3dCSCG_Trace_Element

import warnings
import matplotlib.pyplot as plt


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

        # let's do the plots ____________________________________________
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


class _3dCSCG_Trace_Elements_SELFCHECK(FrozenOnly):
    """We find some specific groups of elements."""

    def __init__(self, trace_elements):
        self._elements_ = trace_elements
        self._freeze_self_()

    def outward_unit_normal_vector(self, *args, **kwargs):
        return self._elements_.___SELFCHECK_outward_unit_normal_vector___(*args, **kwargs)














if __name__ == '__main__':
    # mpiexec -n 12 python _3dCSCG\mesh\trace\elements\main.py
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