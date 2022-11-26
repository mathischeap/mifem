# -*- coding: utf-8 -*-
import random

from root.config.main import *
from components.freeze.main import FrozenClass
from components.quadrature import Quadrature
from components.exceptions import ElementsLayoutError

# noinspection PyUnresolvedReferences
class CSCG_MESH_BASE(FrozenClass):
    """"""


    def __repr__(self):
        """"""
        return f"{self.ndim}dCSCG-mesh:{self.domain.domain_input.domain_name}={self.elements.GLOBAL_num}-elements"

    def ___PRIVATE_BASE_analyze_element_layout___(self, element_layout):
        """Here we return the Element_Layout (EL) for a particular regions.

        We will return EL which will be used as the element_layout of a particular regions.

        Parameters
        ----------
        element_layout : None, int, tuple, list,
            When it is `None`, we return EL = [1,1,1] (for 3d) or [1,1] (for 2d).

            When it is int, for example, element_layout=i, we return EL = [i,i,i] or [i,i].

            When it is tuple, list, we study each elements_layout[i]

                if elements_layout[i] is string, then we need to parse it.

                    it can be:
                        elements_layout[i]='Lobatto:X',
                        we study it according to Lobatto distribution and get the
                        correct `_el_[i]`

                else:
                    `_el_[i] = elements_layout[i]`

        """

        if element_layout is None:
            _el_ = tuple([1 for _ in range(self.ndim)])

        elif isinstance(element_layout, int):
            assert element_layout >= 1, " <Mesh> : needs more than 1 element."
            _el_ = tuple([element_layout for _ in range(self.ndim)])

        elif element_layout.__class__.__name__ in ('list', 'tuple'):
            assert len(element_layout) == self.ndim, " <Mesh> : element_layout error."
            _el_ = [None for _ in range(self.ndim)]
            for i in range(self.ndim):

                if isinstance(element_layout[i], str):

                    # special element_layout input, using str ..........

                    if element_layout[i][:7] == 'Lobatto':
                        degree = int(element_layout[i].split(':')[1])
                        spacing = Quadrature(degree, category='Lobatto').quad[0]
                        _el_[i] = spacing[1:] - spacing[:-1]
                    else:
                        raise Exception('Element_layout key={} wrong.'.format(element_layout[i]))


                else:
                    _el_[i] = element_layout[i]


        else:
            raise ElementsLayoutError()

        return _el_


    def ___PRIVATE_BASE_get_region_elements_distribution_type___(self):
        """REDT: region_elements_distribution_type

        `self._REDT_` is used to indicate how elements (amount-wise) are distributed to regions. So we do not care about
        the topology or whatever, we only care about the amount of elements in each region:

            1). When the amount of elements in each region is more or less the same,
                self._REDT_ = 'equal'
            else).
                self._REDT_ = 'unknown'

        :return: `self._REDT_`

        """
        rns = self.domain.regions.names
        NE_IR = self._num_elements_in_region_
        NT_TE = self._num_total_elements_
        assert sum(list(NE_IR.values())) == NT_TE, "elements regions-wise distribution is wrong."
        region_wise_ratio: Dict[str] = dict()
        for rn in NE_IR: region_wise_ratio[rn] = NE_IR[rn] / NT_TE
        sf = 1 / len(rns)
        # self._region_wise_element_ratio_ = region_wise_ratio
        self._region_relative_size_: Dict[str] = dict()
        FACTOR = list()
        self.___EDF___ = dict()
        for rn in rns:
            factor = region_wise_ratio[rn] / sf
            self._region_relative_size_[rn] = factor
            FACTOR.append(factor)
            self.___EDF___[rn] = factor
        FACTOR = np.array(FACTOR)

        maxF = np.max(FACTOR)
        minF = np.min(FACTOR)

        # -------------------------------------------------------------------------------------------

        if 0.8 < minF and maxF < 1.2:
            self._REDT_ = 'equal' # REDT: region_elements_distribution_type
        else:
            self._REDT_ = 'unknown'  # REDT: region_elements_distribution_type

        if self.domain.regions.num == 1:
            assert self._REDT_ == 'equal', "This must be the case when #regions == 1."




    def ___PRIVATE_BASE_decide_EDM___(self, EDM):
        """
        We parse which EDM (`self._EDM_`) we finally should use.

        When input `EMD` is in `customized_methods`. Then we just use it. Otherwise, we use our scheme to find the most
        proper 'self._EDM_' (I hope so).

        :param EDM:
        :return:
        """
        customized_methods = ('debug', 'SWV0', 'chaotic')

        if EDM in customized_methods: # we strongly use the method.

            if EDM == 'debug':
                EDM = None

            self._EDM_ = EDM

        elif EDM is None: # we try to find a proper one.

            if self._REDT_ == 'unknown':
                self._EDM_ = None # when Region_Elements_Distribution_Type is unknown, we use trivial element division.

            else:
                if SIZE == 1: # only one core, then there is nothing we can do.
                    self._EDM_ = None
                else:
                    if self._REDT_ == 'equal':
                        if SIZE <= self.domain.regions.num: # we always do this when we have many regions!
                            EDM = 'cores_no_more_than_regions'
                        else:
                            EDM = 'SWV0'
                    else:
                        raise NotImplementedError(f"Can not AUTO-decide EDM for "
                                                  f"Region_Elements_Distribution_Type: {self._REDT_}.")

                    self._EDM_ = EDM

        else:
            raise NotImplementedError(f"EDM = {EDM} not coded.")




    def ___PRIVATE_BASE_parse_element_distribution_method___(self):
        """
        With this method, we must define properties:
            _num_local_elements_, _element_distribution_ and _element_indices_.

        :return:
        """
        # ... we find the three attributes through different methods
        EDM = self._EDM_

        master_core_load_factor = 0.75  # change this factor to change the load of master core.
        secretary_core_load_factor = 0.9 # change this factor to change the load of secretary core.
        # master_core_load_factor to be the larger, master core will be more busy.
        assert 0.1 <= master_core_load_factor <= 0.95, "master_core_load_factor need be in [0.1,0.95]"
        # secretary_core_load_factor to be the larger, secretary core will be more busy.
        assert master_core_load_factor < secretary_core_load_factor <=1, \
            f"secretary_core_load_factor need be in [master_core_load_factor({master_core_load_factor}),1]"

        self.___MC_LF___ = master_core_load_factor
        self.___SC_LF___ = secretary_core_load_factor

        nrs = self.domain.regions.num
        rns = self.domain.regions.names

        disDict = dict()  # keys are core numbers, and values are #elements in this core.

        # with this if-else, we have to give the disDict,
        if EDM is None: # default method

            numOfTotalElements = self._num_total_elements_

            if SIZE == 1:
                disDict[0] = numOfTotalElements
            elif SIZE == 2:
                baseNum = numOfTotalElements // 2
                for i in range(SIZE):
                    if i == MASTER_RANK:
                        disDict[i] = baseNum
                    else:
                        disDict[i] = numOfTotalElements - baseNum
            else:
                numEleForMaster = numOfTotalElements // SIZE
                if numEleForMaster == 1:
                    pass
                elif 5 >= numEleForMaster > 2:
                    numEleForMaster = 2
                else:
                    # change this factor to change the payload for master core.
                    numEleForMaster = int(master_core_load_factor * numEleForMaster)
                disDict[MASTER_RANK] = numEleForMaster
                if SECRETARY_RANK != MASTER_RANK:
                    numEleForSecretary = (numOfTotalElements-numEleForMaster) // (SIZE - 1)
                    if numEleForSecretary >= 7:
                        # change this factor to change the payload for secretary core.
                        numEleForSecretary = int(secretary_core_load_factor * numEleForSecretary)
                    disDict[SECRETARY_RANK] = numEleForSecretary
                    num = numOfTotalElements - numEleForMaster - numEleForSecretary
                    parts = SIZE - 2
                else:
                    num = numOfTotalElements - numEleForMaster
                    parts = SIZE - 1
                eleDis = [num // parts + (1 if x < num % parts else 0) for x in range(parts)]
                for i in range(SIZE):
                    if i not in disDict:
                        assert i in SLAVE_RANKS, "Error."
                        disDict[i] = eleDis.pop(0)
                assert eleDis == list(), "Error."

        elif EDM == 'chaotic':  # randomly distribute the elements (for testing purpose.)
            if RANK == MASTER_RANK:

                while 1:
                    DIS = list()
                    for i in range(SIZE):

                        _  = random.randint(0,100)
                        if _ <= 5:
                            _ = 0

                        DIS.append(_)

                    SUM = sum(DIS)
                    if SUM > 0:
                        break


                DIS = np.array(DIS) / SUM

                np.testing.assert_almost_equal(np.sum(DIS), 1)

                numOfTotalElements = self._num_total_elements_

                _dis_ = list()
                for i in range(SIZE):
                    _dis_.append(int(numOfTotalElements * DIS[i]))

                _dis_[-1] = numOfTotalElements - sum(_dis_[:-1])

                assert sum(_dis_) == numOfTotalElements

                for i in range(SIZE):
                    disDict[i] = _dis_[i]

            else:
                pass

            disDict = COMM.bcast(disDict, root=MASTER_RANK)


        elif EDM == 'cores_no_more_than_regions': # few cores many regions cases.
            assert SIZE > 1, "When only have one core, we use EDM = None!"
            assert nrs >= SIZE, "a trivial check!"

            if self._REDT_ == 'equal':
                if nrs % SIZE == 0: # regions can be equally distributed to cores.

                    regions_per_core = int(nrs / SIZE) # each core will handle this many regions.
                    assert regions_per_core * SIZE == nrs

                    for i in range(SIZE):
                        dDi = 0
                        for j in range(regions_per_core):
                            m = j + i * regions_per_core
                            rn = rns[m]
                            dDi += self._num_elements_in_region_[rn]

                        disDict[i] = dDi

                else:

                    region_distribution = [nrs // SIZE + (1 if x < nrs % SIZE else 0) for x in range(SIZE)][::-1]

                    rns = list(self.domain.regions.names)[::-1]

                    core_regions_dict = dict()
                    for i in range(SIZE):
                        core_regions_dict[i] = list()
                        for j in range(region_distribution[i]):
                            core_regions_dict[i].append(rns.pop())

                    for i in range(SIZE):
                        disDict[i] = 0
                        cover_regions = core_regions_dict[i]
                        for cr in cover_regions:
                            disDict[i] += self._num_elements_in_region_[cr]

            else:
                raise NotImplementedError(f"EDM: 'cores_no_more_than_regions' can not handle type: '{self._REDT_}' "
                                          f"regions-wise element distribution.")

        elif EDM == 'SWV0': # smart way version 0; cores do not have elements from different regions.
            assert SIZE > 1, f"When only have one core, we use EDM = None!"
            assert SIZE >= nrs, f"EDM SWV0 only fits when Num of Cores is not lower than Mum of Regions, NCS={SIZE} < {nrs}=NRS."

            self.___region_cores_dict___ = dict()

            if self._REDT_ == 'equal': # all regions have roughly same amount of elements.
                if nrs == 1 or SIZE % nrs == 0: # cores are equally distributed.

                    cores_per_region = int(SIZE / nrs)

                    for i in range(nrs):
                        rn = rns[i]
                        self.___region_cores_dict___[rn] = list()

                        elements_2b_divided = self._num_elements_in_region_[rn]

                        eleDis = [elements_2b_divided // cores_per_region +
                                  (1 if x < elements_2b_divided % cores_per_region else 0)
                                  for x in range(cores_per_region)]

                        for j in range(cores_per_region):
                            m = j + cores_per_region * i
                            disDict[m] = eleDis[j]

                            self.___region_cores_dict___[rn].append(m)

                else:

                    cores_distribution = [SIZE // nrs + (1 if x < SIZE % nrs else 0) for x in range(nrs)]

                    region_cores_dict = dict()
                    core_list = [_ for _ in range(SIZE)][::-1]

                    for i, rn in enumerate(rns):
                        region_cores_dict[rn] = list()
                        for j in range(cores_distribution[i]):
                            region_cores_dict[rn].append(core_list.pop())

                    self.___region_cores_dict___ = region_cores_dict

                    for rn in region_cores_dict:
                        cores = region_cores_dict[rn]
                        n_l_c = len(cores)
                        e_2_d = self._num_elements_in_region_[rn]
                        eleDis = [e_2_d // n_l_c + (1 if x < e_2_d % n_l_c else 0) for x in range(n_l_c)][::-1]
                        for i, c in enumerate(cores):
                            disDict[c] = eleDis[i]

            else:
                raise NotImplementedError(f"EDM: 'SWV0' can not handle type: '{self._REDT_}' "
                                          f"regions-wise element distribution.")


            master_elements = disDict[MASTER_RANK]
            for rn in self.___region_cores_dict___:
                if MASTER_RANK in self.___region_cores_dict___[rn]:

                    if len(self.___region_cores_dict___[rn]) == 1:
                        pass
                    else:
                        slaves = list()
                        for i in self.___region_cores_dict___[rn]:
                            if i != MASTER_RANK:
                                slaves.append(i)

                        if master_elements <= 2:
                            pass
                        else:
                            if master_elements == 3:
                                master_elements = 2
                                give_to_slaves = 1
                            else:
                                _ = int(master_elements * master_core_load_factor)
                                give_to_slaves = master_elements - _
                                master_elements = _

                            parts = len(slaves)
                            _ = [give_to_slaves // parts + (1 if x < give_to_slaves % parts else 0) for x in range(parts)]

                            for i, SLA in enumerate(slaves):
                                disDict[SLA] += _[i]

                            disDict[MASTER_RANK] = master_elements

                else:
                    pass

            secretary_elements = disDict[SECRETARY_RANK]
            for rn in self.___region_cores_dict___:
                if SECRETARY_RANK in self.___region_cores_dict___[rn]:

                    if len(self.___region_cores_dict___[rn]) == 1:
                        pass
                    elif len(self.___region_cores_dict___[rn]) == 2 and MASTER_RANK in self.___region_cores_dict___[rn]:
                        pass
                    else:
                        slaves = list()
                        for i in self.___region_cores_dict___[rn]:
                            if i != MASTER_RANK and i != SECRETARY_RANK:
                                slaves.append(i)

                        if secretary_elements <= 2:
                            pass
                        else:
                            if secretary_elements == 3:
                                secretary_elements = 2
                                give_to_slaves = 1
                            else:
                                _ = int(secretary_elements * secretary_core_load_factor)
                                give_to_slaves = secretary_elements - _
                                secretary_elements = _

                            parts = len(slaves)
                            _ = [give_to_slaves // parts + (1 if x < give_to_slaves % parts else 0) for x in range(parts)]

                            for i, SLA in enumerate(slaves):
                                disDict[SLA] += _[i]

                            disDict[SECRETARY_RANK] = secretary_elements

            assert self.___region_cores_dict___ != dict()
            core_pool = set()
            for rn in self.___region_cores_dict___:
                core_pool.update(self.___region_cores_dict___[rn])
            assert len(core_pool) == SIZE




        #---- IF new EDM added we must add `elif` below to make disDict for this new EDM !!!!!!











        else:
            raise Exception(f"element_distribution_method: '{EDM}' not coded for "
                            f"<parse_element_distribution_method>.")


        # must not define _num_local_elements_, _element_distribution_, _element_indices_ yet
        assert not hasattr(self, '_num_local_elements_')
        assert not hasattr(self, '_element_distribution_')
        assert not hasattr(self, '_element_indices_')


        # do the things: Now, we must have disDict.
        assert isinstance(disDict, dict)
        for i in range(SIZE): assert i in disDict, f"disDict not full, miss distribution for core #{i}."

        # now, we parse disDict to obtain  _num_local_elements_, _element_distribution_, _element_indices_
        self._num_local_elements_ = disDict[RANK]
        before_elements = 0
        self._element_distribution_ = dict()
        for i in range(0, SIZE):
            self._element_distribution_[i] = range(before_elements, before_elements+disDict[i])
            if i == RANK:
                self._element_indices_ = range(before_elements, before_elements+self._num_local_elements_)
            before_elements += disDict[i]




    def ___PRIVATE_BASE_analyze_element_distribution___(self):
        """"""
        # Distributed! need have defined _num_local_elements_, _element_distribution_, _element_indices_
        assert hasattr(self, '_num_local_elements_')
        assert hasattr(self, '_element_distribution_')
        assert hasattr(self, '_element_indices_')
        # Now do some more checks ...
        check_num_elements = 0
        for r in range(SIZE):
            RANGE = self._element_distribution_[r]
            assert RANGE.start <= RANGE.stop, f"element_distribution in core {r} is wrong."
            check_num_elements += len(RANGE)
            if r == RANK:
                assert self._num_local_elements_ == len(RANGE), f"element_distribution in core {RANK} is wrong."
                assert RANGE == self._element_indices_, f"element_indices in core {RANK} is wrong."
        assert check_num_elements == self._num_total_elements_, f"element_distribution is wrong."
        # ......
        assert self._element_indices_ == self._element_distribution_[RANK]


        # ! a check: we make sure that elements are distributed into cores in an increasing sequence (MUST BE). ------ !
        for i in range(SIZE):
            if RANK == i:
                RANGE = self._element_distribution_[RANK]
                start, stop = RANGE.start, RANGE.stop

                if RANK > 0: # not the first core
                    START = COMM.recv(source=RANK - 1, tag=RANK - 1)

                    assert start == START

                if RANK < SIZE-1: # not the last core
                    COMM.send(stop, dest=RANK + 1, tag=RANK)

        if RANK == SIZE - 1 and SIZE > 1:
            # noinspection PyUnboundLocalVariable
            assert stop == self._num_total_elements_
        # ! ... end check here ----------------------------------------------------------


        ___is_occupying_all_cores___ = \
            all([len(self._element_distribution_[c]) for c in self._element_distribution_])

        if not hasattr(self, '___is_occupying_all_cores___'):
            # if all cores have elements to play with?
            self.___is_occupying_all_cores___ = ___is_occupying_all_cores___
        else:
            assert self.___is_occupying_all_cores___ == ___is_occupying_all_cores___, "Can not change this property."




        ___USEFUL_regions_and_boundaries___ = list()
        _elements_in_regions_ = list()
        region_map = self.domain.regions.map

        for i in self._element_indices_:
            rn = self.do.find.region_name_of_element(i)
            if rn not in ___USEFUL_regions_and_boundaries___:
                ___USEFUL_regions_and_boundaries___.append(rn)

            if rn not in _elements_in_regions_:
                _elements_in_regions_.append(rn)

            for side_rn in region_map[rn]:
                if side_rn not in ___USEFUL_regions_and_boundaries___:
                    ___USEFUL_regions_and_boundaries___.append(side_rn)

            if SAFE_MODE:
                assert RANK == self.DO.FIND_slave_of_element(i)

        for rg in self.domain.domain_input.periodic_boundaries_involved_regions:
            # we just put all regions have periodic boundaries into useful regions. Bad but simple.
            if rg not in ___USEFUL_regions_and_boundaries___:
                ___USEFUL_regions_and_boundaries___.append(rg)




        if not hasattr(self, '___USEFUL_regions_and_boundaries___'):
            # this core will need to know element numbering in these regions
            self.___USEFUL_regions_and_boundaries___ = ___USEFUL_regions_and_boundaries___
        else:
            assert self.___USEFUL_regions_and_boundaries___ == ___USEFUL_regions_and_boundaries___, \
                "Cannot change this property"



        if not hasattr(self, '_elements_in_regions_'):
            # this core's elements are in this region (or these regions)
            self._elements_in_regions_ = _elements_in_regions_
        else:
            assert self._elements_in_regions_ == _elements_in_regions_, "Can not change this property"