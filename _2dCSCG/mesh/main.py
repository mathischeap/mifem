


from root.config.main import *
from inheriting.CSCG.mesh.base import CSCG_MESH_BASE
from screws.decorators.accepts import accepts, memoize5
from screws.exceptions import ElementsLayoutError, ElementEdgePairError
from typing import Dict, Union
from _2dCSCG.mesh.elements.main import _2dCSCG_Mesh_Elements
from _2dCSCG.mesh.trace.main import _2dCSCG_Trace
from _2dCSCG.mesh.visualize.main import _2dCSCG_Mesh_Visualize
from _2dCSCG.mesh.boundaries.main import _2dCSCG_Mesh_Boundaries
from _2dCSCG.mesh.periodic_setting.main import _2dCSCG_PeriodicDomainSetting
from _2dCSCG.mesh.deprecated.coordinate_transformation import CoordinateTransformation
from _2dCSCG.mesh.do.main import _2dCSCG_Mesh_DO


class _2dCSCG_Mesh(CSCG_MESH_BASE):
    """The 3dCSCG mesh."""
    def __init__(self, domain, element_layout=None, EDM=None):
        assert domain.ndim == 2, " <Mesh> "
        self._domain_ = domain
        cOmm.barrier() # for safety reasons

        self._DO_ = _2dCSCG_Mesh_DO(self)
        self.___PRIVATE_parse_element_layout___(element_layout)


        self.___PRIVATE_BASE_get_region_elements_distribution_type___()
        self.___PRIVATE_BASE_decide_EDM___(EDM)
        self.___PRIVATE_BASE_parse_element_distribution_method___()
        self.___PRIVATE_BASE_analyze_element_distribution___()

        self.___character_num_elements___ = dict()
        self.___PRIVATE_generate_element_global_numbering___()
        self.___PRIVATE_optimize_element_distribution___()


        self.___PRIVATE_generate_element_map___()
        self.___PRIVATE_modify_elements_map_wr2_periodic_setting___()
        self.___PRIVATE_generate_boundary_element_edges___()


        self.___DEPRECATED_ct___ = CoordinateTransformation(self) # only for test purpose
        self._elements_ = _2dCSCG_Mesh_Elements(self)
        self._trace_ = _2dCSCG_Trace(self)
        self._visualize_ = _2dCSCG_Mesh_Visualize(self)
        self._boundaries_ = _2dCSCG_Mesh_Boundaries(self)
        self.___define_parameters___ = None
        self.___TEST_MODE___ = False
        self.do.reset_cache()
        self._freeze_self_()




    @accepts('self', (tuple, list, int, dict, 'NoneType'))
    def ___PRIVATE_parse_element_layout___(self, element_layout):
        """We parse the element_layout.

        element_layout should be a dict, if it is not, we make it a dict with same values.
        """
        rns = self.domain.regions.names
        # __ prepare the element_layout _______________________________________________
        if not isinstance(element_layout, dict):
            EL = dict()
            for rn in rns:
                EL[rn] = element_layout
        else:
            EL = element_layout

        self._element_layout_ = dict()
        self._element_ratio_ = dict()
        self._element_spacing_ = dict()
        self._num_elements_in_region_: Dict[str] = dict()
        for rn in rns:
            self._element_layout_[rn], self._element_ratio_[rn], \
            self._element_spacing_[rn], self._num_elements_in_region_[rn] = \
                self.___PRIVATE_parse_element_layout_each_region___(EL[rn])

        self._num_elements_accumulation_ = dict()
        self._num_total_elements_ = 0
        for rn in rns:
            self._num_total_elements_ += self._num_elements_in_region_[rn]
            self._num_elements_accumulation_[self._num_total_elements_] = rn

        self._element_local_numbering_ = dict()
        for rn in rns:
            self._element_local_numbering_[rn] = np.array(
                range(self._num_elements_in_region_[rn])).reshape(
                self._element_layout_[rn], order='F')

    def ___PRIVATE_parse_element_layout_each_region___(self, element_layout):
        """"""

        _el_ = self.___PRIVATE_BASE_analyze_element_layout___(element_layout)

        # ____ check `_el_`, `_el_` is the only one which is used further ______________
        assert len(_el_) == 2
        for i in range(2):
            if isinstance(_el_[i], int):
                assert _el_[i] > 0
            elif _el_[i].__class__.__name__ in ('tuple', 'list', 'ndarray'):
                # _el_[i] represent element_ratio[i]
                assert np.ndim(_el_[i]) == 1, " <Mesh> : elements_layout[{}]={} is wrong.".format(i, _el_[i])
                assert np.min(_el_[i]) > 0, " <Mesh> : elements_layout[{}]={} is wrong.".format(i, _el_[i])
            else:
                raise ElementsLayoutError(" <Mesh> : elements_layout[{}]={} is wrong.".format(i, _el_[i]))
        # We then parse _element_layout_, _element_ratio_, _element_spacing_ __________
        _element_layout_ = [None for _ in range(self.ndim)]
        _element_ratio_ = [None for _ in range(self.ndim)]
        _element_spacing_ = [None for _ in range(self.ndim)]
        for i in range(self.ndim):
            if isinstance(_el_[i], int):
                assert _el_[i] >= 1, " <Mesh> : elements_layout[{}]={} is wrong.".format(i, _el_[i])
                _element_layout_[i] = _el_[i]
                _element_ratio_[i] = 1 / _el_[i] * np.ones(_el_[i])
            elif _el_[i].__class__.__name__ in ('tuple', 'list', 'ndarray'):
                # noinspection PyTypeChecker,PyUnresolvedReferences
                _element_layout_[i] = np.size(_el_[i])
                _element_ratio_[i] = np.array(_el_[i]) / np.sum(np.array(_el_[i]))
            else:
                raise ElementsLayoutError(" <Mesh> : elements_layout[{}]={} is wrong.".format(i, _el_[i]))
            # noinspection PyTypeChecker
            _element_spacing_[i] = np.zeros(_element_layout_[i] + 1)
            # noinspection PyTypeChecker,PyUnresolvedReferences
            _element_spacing_[i][-1] = 1
            # noinspection PyTypeChecker
            for j in range(1, _element_layout_[i]):
                # noinspection PyTypeChecker,PyUnresolvedReferences
                _element_spacing_[i][j] = np.sum(_element_ratio_[i][0:j])
        # Now we some properties_______________________________________________________
        _element_layout_ = tuple(_element_layout_)
        _element_ratio_ = tuple(_element_ratio_)
        _element_spacing_ = tuple(_element_spacing_)
        _num_elements_in_region_ = np.prod(_element_layout_)
        return _element_layout_, _element_ratio_, _element_spacing_, _num_elements_in_region_





    def ___PRIVATE_generate_element_global_numbering___(self, number_what=None):
        """"""

        rns = self.domain.regions.names

        if number_what is None:
            DO_number_what = self.___USEFUL_regions_and_boundaries___
        elif number_what == 'all regions':
            DO_number_what = rns
        elif isinstance(number_what, str) and number_what in rns:
            DO_number_what = [number_what,]
        else:
            raise Exception()

        EDM = self._EDM_

        ___element_global_numbering___ = dict()

        if EDM is None:
            current_num = 0
            for rn in rns:
                if rn in DO_number_what:
                    ___element_global_numbering___[rn] = \
                        np.arange(current_num,
                                  current_num + self._num_elements_in_region_[rn]).reshape(
                                    self._element_layout_[rn], order='F')
                else:
                    pass
                current_num += self._num_elements_in_region_[rn]

        elif EDM == "chaotic":
            current_num = 0
            for rn in rns:
                if rn in DO_number_what:

                    ___element_global_numbering___[rn] = \
                        np.arange(current_num,
                                  current_num + self._num_elements_in_region_[rn]).reshape(
                                    self._element_layout_[rn], order='C')

                else:
                    pass
                current_num += self._num_elements_in_region_[rn]

        elif EDM == 'cores_no_more_than_regions':
            current_num = 0
            for rn in rns:
                if rn in DO_number_what:
                    ___element_global_numbering___[rn] = \
                        np.arange(current_num,
                                  current_num + self._num_elements_in_region_[rn]).reshape(
                                    self._element_layout_[rn], order='F')
                else:
                    pass
                current_num += self._num_elements_in_region_[rn]

        elif EDM == 'SWV0': # smart way version 0; cores do not contain elements from different regions.

            self.___SWV0_para___ = dict()  # once this method will need optimization, we initialize a variable like this.

            current_num = 0

            for rn in rns:
                if rn in DO_number_what:

                    cores_for_this_region = self.___region_cores_dict___[rn]

                    NCR = len(cores_for_this_region)

                    if NCR == 1:

                        EGN = np.arange(current_num,
                                        current_num + self._num_elements_in_region_[rn]).reshape(
                                        self._element_layout_[rn], order='F')

                    else:

                        if NCR < 4: # number of core in this region is 2 or 3.
                            # noinspection PyTupleAssignmentBalance
                            I, J = self._element_layout_[rn]
                            A = [I, J]
                            A.sort()
                            _E_ = np.arange(current_num,
                                            current_num + self._num_elements_in_region_[rn]).reshape(A, order='F')

                        else:

                            # in this region, the element numbering will be range(start, end).
                            start = current_num
                            end = current_num + self._num_elements_in_region_[rn]

                            #!-Once we use _element_distribution_, _element_indices_, or _num_local_elements_, wo do this
                            if rn in self.___character_num_elements___:
                                # to make sure that after optimization, the does not change
                                character_num_elements = self.___character_num_elements___[rn]
                            else:
                                character_num_elements = list()

                                for core in cores_for_this_region:
                                    character_num_elements.append(len(self._element_distribution_[core]))

                                character_num_elements = int(np.mean(character_num_elements))

                                if character_num_elements <= 0: character_num_elements = 1

                                assert character_num_elements <= self._num_elements_in_region_[rn]

                                self.___character_num_elements___[rn] = character_num_elements


                            #!------------------------------------------------------------------------------------------!

                            if character_num_elements <= 3:

                                _E_ = np.arange(current_num,
                                                current_num + self._num_elements_in_region_[rn]).reshape(
                                                self._element_layout_[rn], order='F')

                            else: # we end up with a situation we do not know how to do a proper numbering.

                                # noinspection PyTupleAssignmentBalance
                                I, J = self._element_layout_[rn]

                                A = [I, J]
                                A.sort()
                                A0, A1 = A

                                if A1 / A0 >= NCR * 0.75: # on A1, we have a lot more elements, so we block the regions along A1

                                    _E_ = np.arange(current_num,
                                                    current_num + self._num_elements_in_region_[rn]).reshape(A, order='F')

                                else:  # we now define a general numbering rule.

                                    CNE = character_num_elements  # we use this number to decide how to divide the regions.

                                    if A0 * A0 <= CNE: # A0 * A0 is significantly low.
                                        _E_ = np.arange(current_num,
                                                        current_num + self._num_elements_in_region_[rn]).reshape(A, order='F')

                                    else:

                                        _E_ = np.empty([A0, A1], dtype=int)

                                        r = A0 / A1

                                        Y = (CNE / r) ** 0.5
                                        X = r * Y

                                        X= int(X)
                                        X = 1 if X == 0 else X

                                        if CNE >= 4 and A0 >= 2:
                                            if X < 2:
                                                X = 2

                                        X = A0 if X > A0 else X

                                        B = A0 // X # can have B blocks along A0
                                        R = A0 % X  # will have R elements resting A0

                                        N = start

                                        for n in range(B):
                                            if n != B - 1:
                                                PLUS = X * A1
                                                pylon = np.arange(N, N + PLUS).reshape((X, A1), order='F')
                                                N += PLUS
                                                _E_[n * X:(n + 1) * X, :] = pylon
                                            else:
                                                PLUS = (X + R) * A1
                                                pylon = np.arange(N, N + PLUS).reshape((X + R, A1), order='F')
                                                N += PLUS
                                                _E_[n * X:, :] = pylon

                                        assert N == end, "Something is wrong!, check above lines."

                        # A general scheme to transpose _E_ into EGN ...

                        ESP = _E_.shape                 # shape of _E_
                        DSP = self._element_layout_[rn] # designed shape

                        E0, E1 = ESP

                        if (E0, E1) == DSP:
                            EGN = _E_
                        elif (E1, E0) == DSP:
                            EGN = _E_.T
                        else:
                            raise Exception("SHOULD NEVER REACH HERE.")

                    # give EGN to dict: ___element_global_numbering___ if this region is numbered.
                    ___element_global_numbering___[rn] = EGN

                else: # this region is not numbered, lets pass.
                    pass

                current_num += self._num_elements_in_region_[rn]

        else:
            raise Exception(f"element_distribution_method: '{EDM}' not coded for "
                            f"<generate_element_global_numbering>.")

        # return or save to self ...
        if number_what is None:
            self.___element_global_numbering___ = ___element_global_numbering___
        elif number_what == 'all regions':
            return ___element_global_numbering___
        elif isinstance(number_what, str) and number_what in rns:
            return ___element_global_numbering___[number_what]
        else:
            raise Exception()

    def ___PRIVATE_generate_ALL_element_global_numbering___(self):
        # rns = self.domain.regions.names
        # ALL_element_global_numbering = dict()
        # current_num = 0
        # for rn in rns:
        #     ALL_element_global_numbering[rn] = \
        #         np.arange(current_num, current_num + self._num_elements_in_region_[rn]).reshape(
        #         self._element_layout_[rn], order='F')
        #     current_num += self._num_elements_in_region_[rn]
        # return ALL_element_global_numbering
        return self.___PRIVATE_generate_element_global_numbering___(number_what='all regions')

    def ___PRIVATE_element_division_and_numbering_quality___(self):
        """find the quality of element division (to cores) and element (regions-wise global) numbering quality.

        :return: A tuple of 2 outputs:

                1. The overall quality (of the whole mesh across all cores.)
                2. The local quality of this core.
        """
        if sIze == 1: return 1, 1

        INTERNAL = 0
        EXTERNAL = 0
        BOUNDARY = 0

        for i in self.elements:
            for j in range(4):
                W = self.elements.map[i][j]

                if isinstance(W, str):
                    BOUNDARY += 1
                else:
                    if W in self.elements:
                        INTERNAL += 1
                    else:
                        EXTERNAL += 1

        loc_qua = (INTERNAL + BOUNDARY) / (self.elements.num * 4)

        I = cOmm.reduce(INTERNAL, root=mAster_rank, op=MPI.SUM)
        E = cOmm.reduce(EXTERNAL, root=mAster_rank, op=MPI.SUM)
        B = cOmm.reduce(BOUNDARY, root=mAster_rank, op=MPI.SUM)

        if rAnk == mAster_rank:
            ALL_FACES = self.elements.GLOBAL_num * 4
            assert I + E + B == ALL_FACES, "Something is wrong."
            QUALITY = (I + B) / ALL_FACES
        else:
            QUALITY = None

        QUALITY = cOmm.bcast(QUALITY, root=mAster_rank)

        return QUALITY, loc_qua




    def ___PRIVATE_optimize_element_distribution___(self):
        """After generating global element numbering, we can further do an optimization to further reduce the element
        edge shearing between cores. This will adjust a bit the element distribution in cores, but should not adjust
        too much.

        :return:
        """
        if self._EDM_ == 'SWV0':

            JUST_PASS = self.___SWV0_para___ == dict()
            JUST_PASS = cOmm.allreduce(JUST_PASS, op=MPI.LAND)

            if JUST_PASS: return

            #TODO: optimize it (But not very necessary for 2-d mesh at all. So, may be just leave it.)

        else: # no need to optimize
            return

        # has to do another check ...
        self.___PRIVATE_BASE_analyze_element_distribution___()




    def ___PRIVATE_generate_element_map___(self):
        """We now study the domain.regions.map to generate elements.map which will be the key property of a mesh,
        because it actually records the topology of a mesh.
        """
        self.___element_map___: Dict[int, Union[tuple, list]] = dict()
        for i in self._element_indices_:
            region_name, local_indices = self.___PRIVATE_do_find_region_name_and_local_indices_of_element___(i)
            _em_i_ = self.___PRIVATE_fetch_side_element___(region_name, local_indices)
            self.___element_map___[i] = _em_i_
        if len(self._element_indices_) == 0: # to make sure we initialized the memoize cache.
            self.___PRIVATE_do_find_region_name_and_local_indices_of_element___(-1)

    def ___PRIVATE_fetch_side_element___(self, region_name, local_indices):
        """We try to find the global numbering of the elements or boundary attaching to the local element indexed
        "local_indices" in the regions named "region_name".

        Parameters
        ----------
        region_name : str
        local_indices : tuple

        """
        _side_element_ = dict()
        # _N___________________________________________________________________________
        _side_element_['U'] = None
        if local_indices[0] == 0:
            # then this element's Upper edge is on the regions Upper edge.
            what_attached = self.domain.regions.map[region_name][
                self.domain.regions(region_name)._edge_name_to_index_('U')]
            if what_attached in self.domain._boundary_names_:
                _side_element_['U'] = what_attached
            else:  # must be another regions
                try:
                    _side_element_['U'] = self._element_global_numbering_[what_attached][
                        -1, local_indices[1]]
                except IndexError:
                    raise Exception(
                        " Base layout dis-match between regions {} & {}".format(
                            region_name, what_attached))
        else:  # then this element's Upper element is another element in this regions
            _side_element_['U'] = self._element_global_numbering_[region_name][
                local_indices[0] - 1, local_indices[1]]
        # _S___________________________________________________________________________
        _side_element_['D'] = None
        # noinspection PyTypeChecker
        if local_indices[0] == self._element_layout_[region_name][0] - 1:
            # then this element's Down edge is on the regions Down edge.
            what_attached = self.domain.regions.map[region_name][
                self.domain.regions(region_name)._edge_name_to_index_('D')]
            if what_attached in self.domain._boundary_names_:
                _side_element_['D'] = what_attached
            else:  # must be another regions
                try:
                    _side_element_['D'] = self._element_global_numbering_[what_attached][
                        0, local_indices[1]]
                except IndexError:
                    raise Exception(
                        " Base layout dis-match between regions {} & {}".format(
                            region_name, what_attached))
        else:  # then this element's Down element is another element in this regions
            _side_element_['D'] = self._element_global_numbering_[region_name][
                local_indices[0] + 1, local_indices[1]]
        # _W___________________________________________________________________________
        _side_element_['L'] = None
        if local_indices[1] == 0:
            # then this element's Left edge is on the regions Left edge.
            what_attached = self.domain.regions.map[region_name][
                self.domain.regions(region_name)._edge_name_to_index_('L')]
            if what_attached in self.domain._boundary_names_:
                _side_element_['L'] = what_attached
            else:  # must be another regions
                try:
                    _side_element_['L'] = self._element_global_numbering_[what_attached][
                        local_indices[0], -1]
                except IndexError:
                    raise Exception(
                        " Base layout dis-match between regions {} & {}".format(
                            region_name, what_attached))
        else:  # then this element's Left element is another element in this regions
            _side_element_['L'] = self._element_global_numbering_[region_name][
                local_indices[0], local_indices[1] - 1]
        # _E___________________________________________________________________________
        _side_element_['R'] = None
        # noinspection PyTypeChecker
        if local_indices[1] == self._element_layout_[region_name][1] - 1:
            # then this element's Right edge is on the regions Right edge.
            what_attached = self.domain.regions.map[region_name][
                self.domain.regions(region_name)._edge_name_to_index_('R')]
            if what_attached in self.domain._boundary_names_:
                _side_element_['R'] = what_attached
            else:  # must be another regions
                try:
                    _side_element_['R'] = self._element_global_numbering_[what_attached][
                        local_indices[0], 0]
                except IndexError:
                    raise Exception(
                        " Base layout dis-match between regions {} & {}".format(
                            region_name, what_attached))
        else:  # then this element's Right element is another element in this regions
            _side_element_['R'] = self._element_global_numbering_[region_name][
                local_indices[0], local_indices[1] + 1]
        # ------------------------------------------------------------------------------
        _se_ = list()
        for i in range(4):
            side_name = self.domain.regions(region_name)._edge_index_to_name_(i)
            _se_.append(_side_element_[side_name])
        return _se_





    def ___PRIVATE_initializing_periodic_setting___(self):
        """"""
        pBPs = self.domain.domain_input.periodic_boundary_pairs

        self._periodic_setting_ = _2dCSCG_PeriodicDomainSetting(self, pBPs)
        CES = list()
        for KEYi in self.periodic_setting.periodic_region_edge_pairs.keys():
            CES.extend(
                self.periodic_setting.periodic_region_edge_pairs[KEYi].correspondance_of_element_edges)
        return CES

    def ___PRIVATE_modify_elements_map_wr2_periodic_setting___(self):
        """"""
        ___USEFUL_periodicElementEdgePairs___ = self.___PRIVATE_initializing_periodic_setting___()
        self.___useful_periodic_element_edge_pairs___ = list()
        sideIndexDict = {'U': 0, 'D': 1, 'L': 2, 'R': 3}
        for eachPair in ___USEFUL_periodicElementEdgePairs___:
            pairType, elements, sides = self.do.parse_element_edge_pair(eachPair)[:3]
            if pairType == 'regular|regular':
                elementOne, elementTwo = elements
                sideOne, sideTwo = sides
                if elementOne in self.___element_map___:
                    self.___element_map___[elementOne][sideIndexDict[sideOne]] = elementTwo
                if elementTwo in self.___element_map___:
                    self.___element_map___[elementTwo][sideIndexDict[sideTwo]] = elementOne
                if elementOne in self.___element_map___ or elementTwo in self.___element_map___:
                    self.___useful_periodic_element_edge_pairs___.append(eachPair)
            else:
                raise ElementEdgePairError(f"Pair: {pairType} is not understandable.")

        for i in self.___element_map___:
            self.___element_map___[i] = tuple(self.___element_map___[i])

    def ___PRIVATE_generate_boundary_element_edges___(self):
        """"""
        self.___boundary_element_edges___ = dict()
        for bn in self.domain._boundary_names_:
            self.___boundary_element_edges___[bn] = ()

        for i in self._element_indices_:
            region_name, local_indices = self.___PRIVATE_do_find_region_name_and_local_indices_of_element___(i)
            _em_i_ = self.___element_map___[i]
            for j in range(len(_em_i_)):
                if _em_i_[j] in self.domain._boundary_names_:
                    side_name = self.domain.regions(region_name)._edge_index_to_name_(j)
                    self.___boundary_element_edges___[_em_i_[j]] += (str(i) + '-' + side_name,)




    @property
    def ___statistic___(self):
        """This is a reserved property. It will be called from property `statistic` which is an inherited `property` for
        any class that inherit `FrozenClass` class.
        """
        _dict_ = dict()
        _dict_['total element number'] = self._num_total_elements_
        return _dict_

    @property
    def ___parameters___(self):
        """
        This `parameters` is used to compare if meshes are the same. Therefore, the
        `___parameters___` should uniquely identify a mesh. We also use it tor save and restore a mesh.

        So it is mandatory for saving a mesh.
        """
        return self.___define_parameters___

    def __eq__(self, other):
        return self.standard_properties.parameters == other.standard_properties.parameters




    def ___PRIVATE_reset_cache___(self):
        self.trace.___PRIVATE_reset_cache___()
        self.elements.___PRIVATE_reset_cache___()
        self.boundaries.___PRIVATE_reset_cache___()
        self.___element_global_numbering___ = None

    @memoize5 # must use memoize
    def ___PRIVATE_do_find_region_name_and_local_indices_of_element___(self, i):
        """ Find the regions and the local numbering of ith element. """
        if i == -1: return None # to make sure we initialized the memoize cache.
        region_name = None
        for num_elements_accumulation in self._num_elements_accumulation_:
            if i < num_elements_accumulation:
                region_name = self._num_elements_accumulation_[num_elements_accumulation]
                break
        try:
            local_indices = tuple(np.argwhere(self._element_global_numbering_[region_name] == i)[0])
        except TypeError:
            # we have _element_global_numbering_ is None, but still try to use it, must be in TEST MODE
            assert self.___TEST_MODE___, 'This happens only when TEST MODE is ON.'
            if hasattr(self, '___TEST_cache___'):
                local_indices = self.___TEST_cache___[i][1]
            else:
                self._melt_self_()
                self.___TEST_cache___ = dict()
                self._freeze_self_()
                AEGN = self.___PRIVATE_generate_ALL_element_global_numbering___()
                for j in range(self._num_total_elements_):
                    if j not in self._element_indices_:
                        rnj = None
                        for num_elements_accumulation in self._num_elements_accumulation_:
                            if j < num_elements_accumulation:
                                rnj = self._num_elements_accumulation_[num_elements_accumulation]
                                break
                        LIj = tuple(np.argwhere(AEGN[rnj] == j)[0])
                        self.___TEST_cache___[j] = [rnj, LIj]
                return self.___TEST_cache___[i]
        return region_name, local_indices

    @property
    def _element_global_numbering_(self):
        return self.___element_global_numbering___




    @property
    def domain(self):
        return self._domain_

    @property
    def ndim(self):
        return self.domain.ndim

    @property
    def do(self):
        return self._DO_

    @property
    def elements(self):
        return self._elements_

    @property
    def periodic_setting(self):
        return self._periodic_setting_

    @property
    def trace(self):
        return self._trace_

    @property
    def visualize(self):
        return self._visualize_

    @property
    def boundaries(self):
        return self._boundaries_