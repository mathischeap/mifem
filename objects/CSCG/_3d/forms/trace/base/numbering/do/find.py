



from screws.freeze.main import FrozenOnly





class _3dCSCG_Trace_Numbering_DO_FIND(FrozenOnly):
    def __init__(self, DO):
        self._numbering_ = DO._numbering_
        self._0TraceLocalCache_ = dict()
        self._1TraceLocalCache_ = dict()
        self._2TraceLocalCache_ = dict()
        self._freeze_self_()

    def local_dofs_on_element_side(self, side_name):
        """

        :param side_name:
        :return:
        """
        numbering = self._numbering_

        if numbering._tf_.k == 0:
            return self.___PRIVATE_find_0Trace_local_dofs_on_element_side___(side_name)
        if numbering._tf_.k == 1:
            return self.___PRIVATE_find_1Trace_local_dofs_on_element_side___(side_name)
        if numbering._tf_.k == 2:
            return self.___PRIVATE_find_2Trace_local_dofs_on_element_side___(side_name)
        else:
            raise NotImplementedError(f"not coded for {numbering._tf_.k}-trace-form.")


    def ___PRIVATE_find_0Trace_local_dofs_on_element_side___(self, side_name):
        """

        :param side_name:
        :return:
        """
        if side_name in self._0TraceLocalCache_: return self._0TraceLocalCache_[side_name]

        nbc = self._numbering_._tf_.num.basis_components
        num_NS = nbc['N'][0]
        num_WE = nbc['W'][0]
        num_BF = nbc['B'][0]

        if side_name == 'N':
            self._0TraceLocalCache_['N'] = [i for i in range(num_NS)]
        elif side_name == 'S':
            self._0TraceLocalCache_['S'] = [i for i in range(num_NS, 2*num_NS)]
        elif side_name == 'W':
            self._0TraceLocalCache_['W'] = [i for i in range(2*num_NS, 2*num_NS+num_WE)]
        elif side_name == 'E':
            self._0TraceLocalCache_['E'] = [i for i in range(2*num_NS+num_WE, 2*num_NS+2*num_WE)]
        elif side_name == 'B':
            self._0TraceLocalCache_['B'] = [i for i in range(2*num_NS+2*num_WE, 2*num_NS+2*num_WE+num_BF)]
        elif side_name == 'F':
            self._0TraceLocalCache_['F'] = [i for i in range(2*num_NS+2*num_WE+num_BF, 2*num_NS+2*num_WE+2*num_BF)]
        else:
            raise Exception()

        return self._0TraceLocalCache_[side_name]

    def ___PRIVATE_find_1Trace_local_dofs_on_element_side___(self, side_name):
        """

        :param side_name:
        :return:
        """
        if side_name in self._1TraceLocalCache_: return self._1TraceLocalCache_[side_name]

        nbc = self._numbering_._tf_.num.basis_components
        num_NS = nbc['N'][0]
        num_WE = nbc['W'][0]
        num_BF = nbc['B'][0]

        if side_name == 'N':
            self._1TraceLocalCache_['N'] = [i for i in range(num_NS)]
        elif side_name == 'S':
            self._1TraceLocalCache_['S'] = [i for i in range(num_NS, 2*num_NS)]
        elif side_name == 'W':
            self._1TraceLocalCache_['W'] = [i for i in range(2*num_NS, 2*num_NS+num_WE)]
        elif side_name == 'E':
            self._1TraceLocalCache_['E'] = [i for i in range(2*num_NS+num_WE, 2*num_NS+2*num_WE)]
        elif side_name == 'B':
            self._1TraceLocalCache_['B'] = [i for i in range(2*num_NS+2*num_WE, 2*num_NS+2*num_WE+num_BF)]
        elif side_name == 'F':
            self._1TraceLocalCache_['F'] = [i for i in range(2*num_NS+2*num_WE+num_BF, 2*num_NS+2*num_WE+2*num_BF)]
        else:
            raise Exception()

        return self._1TraceLocalCache_[side_name]

    def ___PRIVATE_find_2Trace_local_dofs_on_element_side___(self, side_name):
        """

        :param side_name:
        :return:
        """
        if side_name in self._2TraceLocalCache_: return self._2TraceLocalCache_[side_name]

        nbc = self._numbering_._tf_.num.basis_components
        num_NS = nbc['N'][0]
        num_WE = nbc['W'][0]
        num_BF = nbc['B'][0]

        if side_name == 'N':
            self._2TraceLocalCache_['N'] = [i for i in range(num_NS)]
        elif side_name == 'S':
            self._2TraceLocalCache_['S'] = [i for i in range(num_NS, 2*num_NS)]
        elif side_name == 'W':
            self._2TraceLocalCache_['W'] = [i for i in range(2*num_NS, 2*num_NS+num_WE)]
        elif side_name == 'E':
            self._2TraceLocalCache_['E'] = [i for i in range(2*num_NS+num_WE, 2*num_NS+2*num_WE)]
        elif side_name == 'B':
            self._2TraceLocalCache_['B'] = [i for i in range(2*num_NS+2*num_WE, 2*num_NS+2*num_WE+num_BF)]
        elif side_name == 'F':
            self._2TraceLocalCache_['F'] = [i for i in range(2*num_NS+2*num_WE+num_BF, 2*num_NS+2*num_WE+2*num_BF)]
        else:
            raise Exception()

        return self._2TraceLocalCache_[side_name]