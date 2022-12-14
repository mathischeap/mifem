# -*- coding: utf-8 -*-
"""Our mesh have the following structures:

Mesh main structure:
    Mesh -> Domain -> Regions
               |
               ---> DomainInput

Extension structures:
    Mesh -> TraceMesh -> TraceElements -> TraceElement
         |
         -> EdgeMesh
         |
         -> NodeMesh

Components:
    Geometry
    Elements: Mesh -> Elements -> Element

"""
import matplotlib.pyplot as plt
from typing import Dict, Union
from root.config.main import RANK, MASTER_RANK, COMM, np, MPI, SIZE, SECRETARY_RANK
from objects.CSCG.base.mesh.base import CSCG_MESH_BASE
from components.decorators.all import accepts, memoize5#, memoize2
from components.exceptions import ElementsLayoutError, ElementSidePairError
from components.miscellaneous.timer import break_list_into_parts
from objects.CSCG._3d.mesh.elements.main import _3dCSCG_Mesh_Elements
from objects.CSCG._3d.mesh.periodicSetting.main import _3dCSCG_PeriodicDomainSetting
from objects.CSCG._3d.mesh.legacy.coordinate_transformation.transformer import \
    CoordinateTransformation as ___DCT___
from objects.CSCG._3d.mesh.visualize.main import _3dCSCG_Mesh_Visualize
from objects.CSCG._3d.mesh.boundaries.main import _3dCSCG_Mesh_Boundaries
from objects.CSCG._3d.mesh.subGeometry.main import _3dCSCG_Mesh_SubGeometry
from objects.CSCG._3d.mesh.do.main import _3dCSCG_Mesh_DO
from objects.CSCG._3d.mesh.whether import _3dCSCG_Mesh_Whether

from objects.CSCG._3d.mesh.node.main import _3dCSCG_Node
from objects.CSCG._3d.mesh.edge.main import _3dCSCG_Edge
from objects.CSCG._3d.mesh.trace.main import _3dCSCG_Trace

class _3dCSCG_Mesh(CSCG_MESH_BASE):
    """The 3dCSCG mesh."""
    def __init__(self, domain, element_layout=None, EDM=None):
        """

        :param domain: **This will already be decided by the domain_inputs.**
        :param element_layout:
            It should be a dict whose keys are the regions names. If it is not a dict, when we will make
            a dict whose values (all the same) are `element_layout`, which means we use the same element_layout
            in all regions.

            Now for example,
                element_layout = {'R:R1': EL1,
                                  ......,
                                  }

                If `EL1` is None, we will make it become EL1 = (1,1,1).
                If `EL1` is an int, for example, EL1 = i, we will make it become EL1 = (i,i,i).

                So in general EL1 will become of format EL1=(a,b,c), Each entry, a, b or c, must be one of
                    (1): An positive int
                    (2): A 1-d array whose entries are all positive int or float.
                    (3): A string: special, we will parse this string then.

        :param EDM:
            Element-Distribution-Methods; it can be one of (None, 'debug', 'SWV0', 'chaotic')

                None (default): We will try to find a proper method. If we cannot, we will use a most rigid method.
                'debug': We will force ourselves to use the most rigid method.
                'SWV0': Smart-Way-Version-0; A not so small way.
                'chaotic': A very chaotic one. We better not save it because, when we re-build it, it will
                    randomly generate the element distribution again. So we will actually get different meshes
                    even we call same amount of cores, which sometimes is OKAY, but sometimes is not.

        """
        assert domain.ndim == 3, " <Mesh> "
        self._domain_ = domain
        COMM.barrier() # for safety reason

        self.___chaotic_EGN_cache___ = dict()

        self._DO_ = _3dCSCG_Mesh_DO(self)
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
        self.___PRIVATE_generate_boundary_element_sides___()

        self.___DEPRECATED_ct___ = ___DCT___(self) # only for test purpose
        self._elements_ = _3dCSCG_Mesh_Elements(self)
        self._trace_ = None
        self._edge_ = None
        self._node_ = None
        self._visualize_ = None
        self._boundaries_ = None
        self._sub_geometry_ = None
        self.___define_parameters___ = None
        self.___TEST_MODE___ = False
        self.___element_global_numbering___ = None
        self._whether_ = None
        self._freeze_self_()

    @accepts('self', (tuple, list, int, dict, "NoneType"))
    def ___PRIVATE_parse_element_layout___(self, element_layout):
        rns = self.domain.regions.names
        # __ prepare the element_layout ...
        if not isinstance(element_layout, dict):
            EL = {}
            for rn in rns:
                EL[rn] = element_layout
        else:
            EL = element_layout

        # We first parse the input ...
        self._element_layout_: Dict[str] = dict()
        self._element_ratio_: Dict[str] = dict()
        self._element_spacing_: Dict[str] = dict()
        self._num_elements_in_region_: Dict[str] = dict()
        for rn in rns:
            self._element_layout_[rn], self._element_ratio_[rn], \
            self._element_spacing_[rn], self._num_elements_in_region_[rn] = \
                self.___PRIVATE_parse_element_layout_each_region___(EL[rn])
        self._num_total_elements_ = 0
        self._num_elements_accumulation_: Dict = dict()
        for rn in rns:
            self._num_total_elements_ += self._num_elements_in_region_[rn]
            self._num_elements_accumulation_[self._num_total_elements_] = rn

    def ___PRIVATE_parse_element_layout_each_region___(self, element_layout):
        """
        When we reach here, each entry of `element_layout` can only be
            (1) a positive int
            (2) a 1d array whose entry is positive.

        So element_layout = (a, b, c), a or b or c can only be one of (1) and (2) above.

        """

        _el_ = self.___PRIVATE_BASE_analyze_element_layout___(element_layout)

        # __ check _el_, nothing but _el_ goes beyond ......

        assert len(_el_) == 3

        for i in range(self.ndim):
            if isinstance(_el_[i], int):
                assert _el_[i] > 0
            elif _el_[i].__class__.__name__ in ('tuple', 'list', 'ndarray'):
                assert np.ndim(_el_[i]) == 1, \
                    " <Mesh> : elements_layout[{}]={} is wrong.".format(i, _el_[i])
                assert np.min(_el_[i]) > 0, \
                    " <Mesh> : elements_layout[{}]={} is wrong.".format(i, _el_[i])
            else:
                raise ElementsLayoutError(
                    " <Mesh> : elements_layout[{}]={} is wrong.".format(i, _el_[i]))


        # We then parse _element_layout_, _element_ratio_, _element_spacing_ ----------
        _element_layout_: list = [None for _ in range(self.ndim)]
        _element_ratio_: list = [None for _ in range(self.ndim)]
        _element_spacing_: list = [None for _ in range(self.ndim)]
        for i in range(self.ndim):
            if isinstance(_el_[i], int):
                assert _el_[i] >= 1, \
                    " <Mesh> : elements_layout[{}]={} is wrong.".format(i, _el_[i])
                _element_layout_[i] = _el_[i]
                _element_ratio_[i] = 1 / _el_[i] * np.ones(_el_[i])
            elif _el_[i].__class__.__name__ in ('tuple', 'list', 'ndarray'):
                _element_layout_[i] = np.size(_el_[i])
                _element_ratio_[i] = np.array(_el_[i]) / np.sum(np.array(_el_[i]))
            else:
                raise ElementsLayoutError(
                    " <Mesh> : elements_layout[{}]={} is wrong.".format(i, _el_[i]))
            _element_spacing_[i] = np.zeros(_element_layout_[i] + 1)
            _element_spacing_[i][-1] = 1
            for j in range(1, _element_layout_[i]):
                _element_spacing_[i][j] = np.sum(_element_ratio_[i][0:j])

        # Now we some properties ...
        _element_layout_: tuple = tuple(_element_layout_)
        _element_ratio_: tuple = tuple(_element_ratio_)
        _element_spacing_: tuple = tuple(_element_spacing_)
        _num_elements_in_region_ = np.prod(_element_layout_)

        return _element_layout_, _element_ratio_, _element_spacing_, _num_elements_in_region_

    def ___PRIVATE_generate_element_global_numbering___(self, number_what=None):
        """
        IMPORTANT: we can number in whatever sequence within a regions, but cross-regions, we must number them in the
        sequence : regions.names.

        :param number_what:
        :return:
        """
        EDM = self._EDM_
        rns = self.domain.regions.names

        if number_what is None:
            DO_number_what = self.___USEFUL_regions_and_boundaries___
        elif number_what == 'all regions':
            DO_number_what = rns
        elif isinstance(number_what, str) and number_what in rns:
            DO_number_what = [number_what,]
        else:
            raise Exception()

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

            self.___SWV0_para___ = dict() # once this method will need optimization, we initialize a variable like this.

            current_num = 0

            for rn in rns:
                if rn in DO_number_what:

                    cores_for_this_region = self.___region_cores_dict___[rn]

                    NCR = len(cores_for_this_region)

                    # !-Once we use _element_distribution_, _element_indices_, or _num_local_elements_, wo do this:   !!
                    if not hasattr(self, '___have_empty_core___'):
                        self.___have_empty_core___ = dict()
                    if rn in self.___have_empty_core___:
                        # to make sure that after optimization, character_num_elements does not change
                        have_empty_core = self.___have_empty_core___[rn]
                    else:
                        have_empty_core = False
                        for c_iii in cores_for_this_region:
                            if len(self._element_distribution_[c_iii]) == 0:
                                have_empty_core = True
                                break
                        self.___have_empty_core___[rn] = have_empty_core
                    # !--------------------------------------------------------------------------------------------   !!

                    if NCR == 1 or have_empty_core:

                        EGN = np.arange(current_num,
                                        current_num + self._num_elements_in_region_[rn]).reshape(
                                        self._element_layout_[rn], order='F')

                    else:

                        if NCR < 4: # number of core in this regions is 2 or 3.
                            I, J, K = self._element_layout_[rn]
                            A = [I, J, K]
                            A.sort()
                            _E_ = np.arange(current_num,
                                            current_num + self._num_elements_in_region_[rn]).reshape(A, order='F')

                        else:

                            # in this regions, the element numbering will be range(start, end).
                            start = current_num
                            end = current_num + self._num_elements_in_region_[rn]

                            #!-Once we use _element_distribution_, _element_indices_, or _num_local_elements_, wo do this
                            if rn in self.___character_num_elements___:
                                # to make sure that after optimization, character_num_elements does not change
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


                            I, J, K = self._element_layout_[rn] # to determine which scheme to do the numbering.

                            if character_num_elements <= 3:

                                _E_ = np.arange(current_num,
                                                current_num + self._num_elements_in_region_[rn]).reshape(
                                                self._element_layout_[rn], order='F')

                            elif I == J == K:
                                # prepare the memory.
                                _E_ = np.empty(self._element_layout_[rn], dtype=int)

                                # same amount elements along all directions
                                for i in range(1, I+2):
                                    if i**3 > character_num_elements:
                                        break
                                # noinspection PyUnboundLocalVariable
                                i -= 1
                                if i == 0: i = 1

                                if character_num_elements >= 4 and i < 2:
                                    i = 2

                                if i > I: i = I
                                # we will number in small block of i*i elements at x-y plain. And go along z-direction.

                                B = I // i # can have B blocks along x and y.
                                R = I % i  # will have R elements resting along x and y.

                                N = start

                                self.___SWV0_para___[rn] = [I,]
                                for n in range(B):
                                    if n != B-1:
                                        for m in range(B):
                                            if m != B-1:
                                                PLUS = i * i * I
                                                pylon = np.arange(N, N + PLUS).reshape((i, i, I), order='F')
                                                N += PLUS
                                                _E_[m*i:(m+1)*i, n*i:(n+1)*i, : ] = pylon
                                                self.___SWV0_para___[rn].append(PLUS)

                                            else:
                                                PLUS = (i+R) * i * I
                                                pylon = np.arange(N, N + PLUS).reshape((i+R, i, I), order='F')
                                                N += PLUS
                                                _E_[m*i:, n*i:(n+1)*i, : ] = pylon
                                                self.___SWV0_para___[rn].append(PLUS)

                                    else:
                                        for m in range(B):
                                            if m != B-1:
                                                PLUS = i * (i+R) * I
                                                pylon = np.arange(N, N + PLUS).reshape((i, i+R, I), order='F')
                                                N += PLUS
                                                _E_[m*i:(m+1)*i, n*i:, : ] = pylon
                                                self.___SWV0_para___[rn].append(PLUS)

                                            else:
                                                PLUS = (i+R) * (i+R) * I
                                                pylon = np.arange(N, N + PLUS).reshape((i+R, i+R, I), order='F')
                                                N += PLUS
                                                _E_[m*i:, n*i:, : ] = pylon
                                                self.___SWV0_para___[rn].append(PLUS)

                                assert N == end, "must be like this!"

                            else: # we end up with a situation we do not know how to do a proper numbering.

                                A = [I, J, K]
                                A.sort()
                                A0, A1, A2 = A

                                if A2 / A1 >= NCR * 0.75: # on A2, we have a lot more elements, so we block the regions along A2

                                    _E_ = np.arange(current_num,
                                                    current_num + self._num_elements_in_region_[rn]).reshape(A, order='F')

                                else: # we now define a general numbering rule.

                                    CNE = character_num_elements  # we use this number to decide how to divide the regions.

                                    if A0 * A1 * A0 <= CNE: # A0 * A1 is significantly low.
                                        _E_ = np.arange(current_num,
                                                        current_num + self._num_elements_in_region_[rn]).reshape(A, order='F')
                                    else:

                                        _E_ = np.empty([A0, A1, A2], dtype=int)

                                        R01 = A0 / A1
                                        R02 = A0 / A2

                                        Y = (R02 * CNE / R01**2)**(1/3)
                                        X = Y * R01

                                        X, Y = int(X), int(Y)
                                        X = 1 if X == 0 else X
                                        Y = 1 if Y == 0 else Y

                                        if CNE > 4 and A0 >= 2 and A1 >= 2:
                                            if X < 2: X = 2
                                            if Y < 2: Y = 2

                                        X = A0 if X > A0 else X
                                        Y = A1 if Y > A1 else Y

                                        B0 = A0 // X # can have B0 blocks along A0
                                        R0 = A0 % X  # will have R0 elements resting A0
                                        B1 = A1 // Y # can have B1 blocks along A1
                                        R1 = A1 % Y  # will have R1 elements resting A1

                                        N = start

                                        self.___SWV0_para___[rn] = [A2,]
                                        for n in range(B1):
                                            if n != B1 - 1:
                                                for m in range(B0):
                                                    if m != B0 - 1:
                                                        PLUS = X * Y * A2
                                                        pylon = np.arange(N, N + PLUS).reshape((X, Y, A2), order='F')
                                                        N += PLUS
                                                        _E_[m * X:(m + 1) * X, n * Y:(n + 1) * Y, :] = pylon
                                                        self.___SWV0_para___[rn].append(PLUS)

                                                    else:
                                                        PLUS = (X + R0) * Y * A2
                                                        pylon = np.arange(N, N + PLUS).reshape((X + R0, Y, A2), order='F')
                                                        N += PLUS
                                                        _E_[m * X:, n * Y:(n + 1) * Y, :] = pylon
                                                        self.___SWV0_para___[rn].append(PLUS)

                                            else:
                                                for m in range(B0):
                                                    if m != B0 - 1:
                                                        PLUS = X * (Y + R1) * A2
                                                        pylon = np.arange(N, N + PLUS).reshape((X, Y + R1, A2), order='F')
                                                        N += PLUS
                                                        _E_[m * X:(m + 1) * X, n * Y:, :] = pylon
                                                        self.___SWV0_para___[rn].append(PLUS)

                                                    else:
                                                        PLUS = (X + R0) * (Y + R1) * A2
                                                        pylon = np.arange(N, N + PLUS).reshape((X + R0, Y + R1, A2), order='F')
                                                        N += PLUS
                                                        _E_[m * X:, n * Y:, :] = pylon
                                                        self.___SWV0_para___[rn].append(PLUS)

                                        assert N == end, "Something is wrong!, check above lines."

                        # A general scheme to transpose _E_ into EGN ...

                        ESP = _E_.shape                 # shape of _E_
                        DSP = self._element_layout_[rn] # designed shape

                        E0, E1, E2 = ESP

                        if (E0, E1, E2) == DSP:
                            EGN = _E_.transpose((0, 1, 2))
                        elif (E0, E2, E1) == DSP:
                            EGN = _E_.transpose((0, 2, 1))
                        elif (E1, E0, E2) == DSP:
                            EGN = _E_.transpose((1, 0, 2))
                        elif (E1, E2, E0) == DSP:
                            EGN = _E_.transpose((1, 2, 0))
                        elif (E2, E0, E1) == DSP:
                            EGN = _E_.transpose((2, 0, 1))
                        elif (E2, E1, E0) == DSP:
                            EGN = _E_.transpose((2, 1, 0))
                        else:
                            raise Exception("SHOULD NEVER REACH HERE.")

                    # give EGN to dict: ___element_global_numbering___ if this region is numbered.
                    ___element_global_numbering___[rn] = EGN

                else: # this regions is not numbered, lets pass.
                    pass

                current_num += self._num_elements_in_region_[rn]

        else:
            raise Exception(f"element_distribution_method: '{EDM}' not coded for "
                            f"<generate_element_global_numbering>.")

        # check element global numbering ......
        current_num = 0
        for rn in rns:

            if rn in ___element_global_numbering___:
                egn_rn = ___element_global_numbering___[rn]
                assert egn_rn.__class__.__name__ == 'ndarray' and np.ndim(egn_rn) == 3, "must be a 3-d array."
                assert np.min(egn_rn) == current_num and \
                       np.max(egn_rn) == current_num + self._num_elements_in_region_[rn] - 1, \
                       f'Element numbering range in regions {rn} is wrong. Cross regions, the overall numbering must be increasing.'
                # this means within a regions, the element numbering can be anything, but overall, it has to be increasing through rns.

                A = np.shape(egn_rn)
                assert A == self._element_layout_[rn], f"___element_global_numbering___[{rn}] shape wrong!"
                assert self._num_elements_in_region_[rn] == np.prod(A), "A trivial check."
                if self._num_elements_in_region_[rn] < 1000:
                    A = egn_rn.ravel('F')
                    A = set(A)
                    assert len(A) == self._num_elements_in_region_[rn]

            current_num += self._num_elements_in_region_[rn]

        # return or save to self ...
        if number_what is None:
            self.___element_global_numbering___ = ___element_global_numbering___
        elif number_what == 'all regions':
            return ___element_global_numbering___
        elif isinstance(number_what, str) and number_what in rns:
            return ___element_global_numbering___[number_what]
        else:
            raise Exception()

    def ___PRIVATE_generate_element_global_numbering_for_region___(self, region_name):
        """generate element numbering for one regions."""
        return self.___PRIVATE_generate_element_global_numbering___(number_what=region_name)

    def ___PRIVATE_generate_ALL_element_global_numbering___(self):
        return self.___PRIVATE_generate_element_global_numbering___(number_what='all regions')

    def ___PRIVATE_element_division_and_numbering_quality___(self):
        """find the quality of element division (to cores) and element (regions-wise global) numbering quality.

        :return: A tuple of 2 outputs:

                1. The overall quality (of the whole mesh across all cores.) 1 is best, 0 is worst.
                2. The local quality of this core.
        """
        if SIZE == 1: return 1, 1

        INTERNAL = 0
        EXTERNAL = 0
        BOUNDARY = 0

        for i in self.elements:
            for j in range(6):
                W = self.elements.map[i][j]

                if isinstance(W, str):
                    BOUNDARY += 1
                else:
                    if W in self.elements:
                        INTERNAL += 1
                    else:
                        EXTERNAL += 1

        loc_qua = (INTERNAL + BOUNDARY) / (self.elements.num * 6)

        I = COMM.reduce(INTERNAL, root=MASTER_RANK, op=MPI.SUM)
        E = COMM.reduce(EXTERNAL, root=MASTER_RANK, op=MPI.SUM)
        B = COMM.reduce(BOUNDARY, root=MASTER_RANK, op=MPI.SUM)

        if RANK == MASTER_RANK:
            ALL_FACES = self.elements.global_num * 6
            assert I + E + B == ALL_FACES, "Something is wrong."
            QUALITY = (I + B) / ALL_FACES
        else:
            QUALITY = None

        QUALITY = COMM.bcast(QUALITY, root=MASTER_RANK)

        return QUALITY, loc_qua

    def ___PRIVATE_matplot_local_elements___(self):
        """

        :return:
        """
        local_elements = dict() # keys are regions name, values are indices of local elements
        for i in self.elements:
            rn, lid = self.___PRIVATE_do_find_region_name_and_local_indices_of_element___(i)
            if rn not in local_elements:
                local_elements[rn] = (list(), list(), list())
            indices = local_elements[rn]
            indices[0].append(lid[0])
            indices[1].append(lid[1])
            indices[2].append(lid[2])

        LOCAL_ELEMENTS = COMM.gather(local_elements, root=MASTER_RANK)
        if RANK == MASTER_RANK:

            for i, LE in enumerate(LOCAL_ELEMENTS):

                num_regions = len(LE)

                dis = '1' + str(num_regions)

                fig = plt.figure(figsize=(6*num_regions, 6))

                for f, rn in enumerate(LE):
                    plot_num = int(dis + str(f+1))

                    ax = fig.add_subplot(plot_num, projection='3d')
                    indices = LE[rn]

                    ax.scatter(*indices, marker='s')
                    ax.set_xlim3d(0, self._element_layout_[rn][0] - 1)
                    ax.set_ylim3d(0, self._element_layout_[rn][1] - 1)
                    ax.set_zlim3d(0, self._element_layout_[rn][2] - 1)
                    ax.set_xlabel("axis-0")
                    ax.set_ylabel("axis-1")
                    ax.set_zlabel("axis-2")
                    ax.set_title(rn)

                plt.suptitle(f"core #{i}")

                plt.show()
                plt.close()

    def ___PRIVATE_optimize_element_distribution___(self):
        """After generating global element numbering, we can further do a optimization to further reduce the element
        side shearing between cores. This will adjust a bit the element distribution in cores, but should not adjust
        too much.

        :return:
        """
        if self._EDM_ == 'SWV0':

            JUST_PASS = self.___SWV0_para___ == dict()
            JUST_PASS = COMM.allreduce(JUST_PASS, op=MPI.LAND)

            if JUST_PASS: return

            # we first merge all ___SWV0_para___ to master ...
            _PA_ = COMM.gather(self.___SWV0_para___, root=MASTER_RANK)
            if RANK == MASTER_RANK:
                PARA = dict()
                for P in _PA_:
                    for pr in P:
                        if pr in PARA:
                            assert P[pr] == PARA[pr], "___SWV0_para___ in different cores must be the same."
                        else:
                            PARA[pr] = P[pr]
                # Now, all ___SWV0_para___ are in PARA ...

                NEW_DIS = dict()

                for rn in self.___region_cores_dict___:

                    if rn in PARA:
                        Loc_Cor = self.___region_cores_dict___[rn]
                        for c in Loc_Cor: assert c not in NEW_DIS, "Safety checker."
                        loc_Par = PARA[rn]
                        layers = loc_Par[0]
                        blocks = loc_Par[1:]

                        num_blocks = len(blocks)
                        num_LocCor = len(Loc_Cor)

                        if num_blocks < 2:
                            pass
                        elif num_LocCor < num_blocks:
                            pass
                        elif layers == 1:
                            pass
                        else:
                            if num_LocCor == num_blocks:
                                for i, c in enumerate(Loc_Cor):
                                    NEW_DIS[c] = blocks[i]
                            else:
                                _d1_ = [num_LocCor // num_blocks + (1 if x < num_LocCor % num_blocks else 0) for x in range(num_blocks)]

                                ___DO___ = True

                                for i, B in enumerate(blocks):

                                    if _d1_[i] > layers:
                                        ___DO___ = False
                                        break

                                if ___DO___:

                                    _d2_ = break_list_into_parts(Loc_Cor, _d1_)

                                    for i, B in enumerate(blocks):

                                        _d3_ = [layers // _d1_[i] + (1 if x < layers % _d1_[i] else 0) for x in range(_d1_[i])][::-1]

                                        num_ele_per_layer = int(B / layers)
                                        assert num_ele_per_layer * layers == B, "Something is wrong."

                                        if MASTER_RANK in _d2_[i]:

                                            OC = list()

                                            for c, L in zip(_d2_[i], _d3_):
                                                if c == MASTER_RANK:
                                                    if L <= 2:
                                                        ML = L
                                                    else:
                                                        ML = int(self.___MC_LF___ * L)

                                                    TO_OTHER = L - ML

                                                else:
                                                    OC.append(c)

                                            # noinspection PyUnboundLocalVariable
                                            if TO_OTHER == 0 or len(OC) == 0:
                                                pass
                                            else:
                                                NOC = len(OC)
                                                DIS_D = [TO_OTHER // NOC + (1 if x < TO_OTHER % NOC else 0) for x in range(NOC)]

                                                for j, c in enumerate(_d2_[i]):
                                                    if c == MASTER_RANK:
                                                        # noinspection PyUnboundLocalVariable
                                                        _d3_[j] = ML
                                                    else:
                                                        _d3_[j] += DIS_D.pop()

                                        if SECRETARY_RANK in _d2_[i]:
                                            OC = list()
                                            for c, L in zip(_d2_[i], _d3_):
                                                if c == MASTER_RANK:
                                                    pass

                                                elif c == SECRETARY_RANK:
                                                    if L <= 2:
                                                        ML = L
                                                    else:
                                                        ML = int(self.___SC_LF___ * L)

                                                    TO_OTHER = L - ML

                                                else:
                                                    OC.append(c)

                                            if TO_OTHER == 0 or len(OC) == 0:
                                                pass
                                            else:
                                                NOC = len(OC)
                                                DIS_D = [TO_OTHER // NOC + (1 if x < TO_OTHER % NOC else 0) for x in range(NOC)]

                                                for j, c in enumerate(_d2_[i]):
                                                    if c == SECRETARY_RANK:
                                                        _d3_[j] = ML
                                                    elif c == MASTER_RANK:
                                                        pass
                                                    else:
                                                        _d3_[j] += DIS_D.pop()

                                        for c, L in zip(_d2_[i], _d3_):

                                            NEW_DIS[c] = L * num_ele_per_layer

                for c in self._element_distribution_:
                    if c not in NEW_DIS:
                        NEW_DIS[c] = len(self._element_distribution_[c])

                # do the things: Now, we must have disDict ...
                assert isinstance(NEW_DIS, dict)
                for i in range(SIZE): assert i in NEW_DIS, f"NEW_DIS not full, miss distribution for core #{i}."
                ED = dict()
                before_elements = 0
                for i in range(SIZE):
                    ED[i] = range(before_elements, before_elements + NEW_DIS[i])
                    before_elements += NEW_DIS[i]

            else:
                ED = None

            ED = COMM.bcast(ED, root=MASTER_RANK)
            self._element_distribution_ = ED
            self._element_indices_ = self._element_distribution_[RANK]
            self._num_local_elements_ = len(self._element_indices_)

        else: # no need to optimize ...
            return

        # has to do another check ...
        self.___PRIVATE_BASE_analyze_element_distribution___()

    def ___PRIVATE_fetch_side_element___(self, region_name, local_indices):
        """
        We try to find the global numbering of the elements or boundary
        attaching to the local element "local_indices" in regions
        "region_name".

        Parameters
        ----------
        region_name : str
        local_indices : tuple

        """
        _side_element_ = dict()
        # _N ...
        _side_element_['N'] = None
        if local_indices[0] == 0:
            # then this element's North side is on the regions North side.
            what_attached = self.domain.regions.map[region_name][
                self.domain.regions(region_name)._side_name_to_index_('N')]
            if what_attached in self.domain._boundary_names_:
                _side_element_['N'] = what_attached
            else:  # must be another regions
                _side_element_['N'] = self._element_global_numbering_[what_attached][
                    -1, local_indices[1], local_indices[2]]
        else:  # then this element's North element another element in this regions
            _side_element_['N'] = self._element_global_numbering_[region_name][
                local_indices[0] - 1, local_indices[1], local_indices[2]]

        # _S ...
        _side_element_['S'] = None
        if local_indices[0] == int(self._element_layout_[region_name][0] - 1):
            # then this element's South side is on the regions South side.
            what_attached = self.domain.regions.map[region_name][
                self.domain.regions(region_name)._side_name_to_index_('S')]
            if what_attached in self.domain._boundary_names_:
                _side_element_['S'] = what_attached
            else:  # must be another regions

                _side_element_['S'] = self._element_global_numbering_[what_attached][
                    0, local_indices[1], local_indices[2]]


        else:  # then this element's South element is another element in this regions
            try:
                _side_element_['S'] = self._element_global_numbering_[region_name][
                    local_indices[0] + 1, local_indices[1], local_indices[2]]
            except IndexError:
                # seems to be the problem in memoize1, using memoize 5 or 2 is likely OK.
                em = f"global elements numbering shape: {self._element_global_numbering_[region_name].shape}, " \
                     f"request indices: {local_indices[0] + 1, local_indices[1], local_indices[2]}. " \
                     f"If this error happens again, back here."
                raise Exception(em)

        # _W ...
        _side_element_['W'] = None
        if local_indices[1] == 0:
            # then this element's West side is on the regions Left side.
            what_attached = self.domain.regions.map[region_name][
                self.domain.regions(region_name)._side_name_to_index_('W')]
            if what_attached in self.domain._boundary_names_:
                _side_element_['W'] = what_attached
            else:  # must be another regions
                _side_element_['W'] = self._element_global_numbering_[what_attached][
                    local_indices[0], -1, local_indices[2]]
        else:  # then this element's West element another element in this regions
            _side_element_['W'] = self._element_global_numbering_[region_name][
                local_indices[0], local_indices[1] - 1, local_indices[2]]

        # _E ...
        _side_element_['E'] = None
        if local_indices[1] == self._element_layout_[region_name][1] - 1:
            # then this element's East side is on the regions Right side.
            what_attached = self.domain.regions.map[region_name][
                self.domain.regions(region_name)._side_name_to_index_('E')]
            if what_attached in self.domain._boundary_names_:
                _side_element_['E'] = what_attached
            else:  # must be another regions
                _side_element_['E'] = self._element_global_numbering_[what_attached][
                    local_indices[0], 0, local_indices[2]]
        else:  # then this element's East element another element in this regions
            _side_element_['E'] = self._element_global_numbering_[region_name][
                local_indices[0], local_indices[1] + 1, local_indices[2]]

        # _B ...
        _side_element_['B'] = None
        if local_indices[2] == 0:
            # then this element's Back side is on the regions Left side.
            what_attached = self.domain.regions.map[region_name][
                self.domain.regions(region_name)._side_name_to_index_('B')]
            if what_attached in self.domain._boundary_names_:
                _side_element_['B'] = what_attached
            else:  # must be another regions
                _side_element_['B'] = self._element_global_numbering_[what_attached][
                    local_indices[0], local_indices[1], -1]
        else:  # then this element's Back element another element in this regions
            _side_element_['B'] = self._element_global_numbering_[region_name][
                local_indices[0], local_indices[1], local_indices[2] - 1]

        # _E ...
        _side_element_['F'] = None
        if local_indices[2] == self._element_layout_[region_name][2] - 1:
            # then this element's Front side is on the regions Right side.
            what_attached = self.domain.regions.map[region_name][
                self.domain.regions(region_name)._side_name_to_index_('F')]
            if what_attached in self.domain._boundary_names_:
                _side_element_['F'] = what_attached
            else:  # must be another regions
                _side_element_['F'] = self._element_global_numbering_[what_attached][
                    local_indices[0], local_indices[1], 0]
        else:  # then this element's Front element another element in this regions
            _side_element_['F'] = self._element_global_numbering_[region_name][
                local_indices[0], local_indices[1], local_indices[2] + 1]

        # ...
        _se_ = list()
        for i in range(6):
            side_name = 'NSWEBF'[i]
            _se_.append(_side_element_[side_name])
        return _se_

    def ___PRIVATE_generate_element_map___(self):
        """We now by studying the self.domain.region_map generate element_map which will be the key property of a mesh,
        because it actually records the topology of a mesh.
        """
        self.___element_map___: Dict[int, Union[tuple, list]] = dict()

        RP = []

        for i in self._element_indices_:
            region_name, local_indices = self.___PRIVATE_do_find_region_name_and_local_indices_of_element___(i)
            _em_i_ = self.___PRIVATE_fetch_side_element___(region_name, local_indices)
            self.___element_map___[i] = _em_i_

            if region_name not in RP:
                RP.append(RP)

        self.___involved_regions___ = RP


        if len(self._element_indices_) == 0: # to make sure we initialized the memoize cache.
            self.___PRIVATE_do_find_region_name_and_local_indices_of_element___(-1)

    def ___PRIVATE_initializing_periodic_setting___(self):

        pBPs = self.domain.domain_input.periodic_boundary_pairs

        self._periodic_setting_ = _3dCSCG_PeriodicDomainSetting(self, pBPs)
        CES = list()
        for KEYi in self.periodic_setting.periodic_region_side_pairs.keys():
            CES.extend(self.periodic_setting.periodic_region_side_pairs[KEYi].correspondence_of_element_sides)
        return CES

    def ___PRIVATE_modify_elements_map_wr2_periodic_setting___(self):
        """"""
        ___USEFUL_periodicElementSidePairs___ = self.___PRIVATE_initializing_periodic_setting___()
        self.___useful_periodic_element_side_pairs___ = list() # will be used for example when generating trace elements

        self.___local_periodic_element_sides___ = list()
        self.___local_periodic_elements___ = list()

        sideIndexDict = {'N': 0, 'S': 1, 'W': 2, 'E': 3, 'B': 4, 'F': 5}
        for eachPair in ___USEFUL_periodicElementSidePairs___:
            pairType, elements, sides = self.___PRIVATE_do_parse_element_side_pair___(eachPair)[:3]
            if pairType == 'regular|regular':
                elementOne, elementTwo = elements
                sideOne, sideTwo = sides
                if elementOne in self.___element_map___:
                    self.___element_map___[elementOne][sideIndexDict[sideOne]] = elementTwo
                if elementTwo in self.___element_map___:
                    self.___element_map___[elementTwo][sideIndexDict[sideTwo]] = elementOne

                if elementOne in self.___element_map___ or elementTwo in self.___element_map___:
                    self.___useful_periodic_element_side_pairs___.append(eachPair)

                if elementOne in self.___element_map___:
                    self.___local_periodic_element_sides___.append(str(elementOne)+sideOne)
                    if elementOne not in self.___local_periodic_elements___:
                        self.___local_periodic_elements___.append(elementOne)
                if elementTwo in self.___element_map___:
                    self.___local_periodic_element_sides___.append(str(elementTwo)+sideTwo)
                    if elementTwo not in self.___local_periodic_elements___:
                        self.___local_periodic_elements___.append(elementTwo)

            else:
                raise ElementSidePairError(f"Pair: {pairType} is not understandable.")

        for i in self.___element_map___:
            self.___element_map___[i] = tuple(self.___element_map___[i])

    def ___PRIVATE_generate_boundary_element_sides___(self):
        self.___boundary_element_sides___ = dict()
        for bn in self.domain._boundary_names_:
            self.___boundary_element_sides___[bn] = ()

        for i in self._element_indices_:
            region_name, local_indices = self.___PRIVATE_do_find_region_name_and_local_indices_of_element___(i)
            _em_i_ = self.___element_map___[i]
            for j in range(len(_em_i_)):
                if _em_i_[j] in self.domain._boundary_names_:
                    side_name = self.domain.regions(region_name)._side_index_to_name_(j)
                    self.___boundary_element_sides___[_em_i_[j]] += (str(i) + '-' + side_name,)

    @property
    def ___statistic___(self):
        """
        This is a reserved property. It will be called from property `statistic` which
        is an inherited `property` for any class that inherit `FrozenClass` class.
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
        if self is other:
            return True
        else:
            return self.standard_properties.parameters == other.standard_properties.parameters

    @staticmethod
    def ___PRIVATE_do_parse_element_side_pair___(eP: str):
        """Element side pairs are also used for trace element keys."""
        if eP.count('-') == 1:
            # must be regular pair to domain boundary!
            elementOne = int(eP.split('-')[0])
            sideOne = eP.split('-')[1][0]
            boundaryName = eP.split('|')[1]
            return 'regular|domainBoundary', elementOne, sideOne, boundaryName
        elif eP.count('-') == 2:
            elementOne, pairTypeINFO, elementTwo = eP.split('-')
            elementOne = int(elementOne)
            elementTwo = int(elementTwo)
            if len(pairTypeINFO) == 3 and pairTypeINFO[1] == '|':
                # regular pair; conforming pair; N|S, W|E, B|F and no twist!
                sideOne = pairTypeINFO[0]
                sideTwo = pairTypeINFO[2]
                return 'regular|regular', \
                       [elementOne, elementTwo], sideOne + sideTwo, None  # None for future extension.
                # for all kinds of pair, return has follow this same rule!
            else:
                raise ElementSidePairError(f"Pair: {pairTypeINFO} is not understandable.")
        else:
            raise Exception('elementSidePair format wrong!')

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
    def do(self):
        return self._DO_

    @property
    def elements(self):
        return self._elements_

    @property
    def periodic_setting(self):
        return self._periodic_setting_

    @property
    def domain(self):
        return self._domain_

    @property
    def ndim(self):
        return self.domain.ndim

    @property
    def trace(self):
        """The trace (face) mesh"""
        if self._trace_ is None:
            self._trace_ = _3dCSCG_Trace(self)
        return self._trace_

    @property
    def edge(self):
        """The edge mesh!"""
        if self._edge_ is None:
            self._edge_ = _3dCSCG_Edge(self)
        return self._edge_

    @property
    def node(self):
        """The node mesh!"""
        if self._node_ is None:
            self._node_ = _3dCSCG_Node(self)
        return self._node_

    @property
    def visualize(self):
        if self._visualize_ is None:
            self._visualize_ = _3dCSCG_Mesh_Visualize(self)
        return self._visualize_

    @property
    def boundaries(self):
        if self._boundaries_ is None:
            self._boundaries_ = _3dCSCG_Mesh_Boundaries(self)
        return self._boundaries_

    @property
    def sub_geometry(self):
        if self._sub_geometry_ is None:
            self._sub_geometry_ = _3dCSCG_Mesh_SubGeometry(self)
        return self._sub_geometry_

    @property
    def quality(self):
        """A factor in [0,1] that reflects the quality of the mesh; 1
        the best, 0 the worst.
        """
        return self.trace.quality['average quality']

    @property
    def whether(self):
        if self._whether_ is None:
            self._whether_ = _3dCSCG_Mesh_Whether(self)
        return self._whether_