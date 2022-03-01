


from screws.frozen import FrozenOnly





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




