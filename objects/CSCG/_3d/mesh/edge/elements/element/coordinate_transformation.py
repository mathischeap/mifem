# -*- coding: utf-8 -*-


import sys
if './' not in sys.path: sys.path.append('./')

from screws.freeze.main import FrozenOnly



class _3dCSCG_Edge_Element_CT(FrozenOnly):
    """"""
    def __init__(self, ee):
        """"""
        self._ee_ = ee
        self._freeze_self_()

    def mapping(self, *ep3, from_element=None, corner_edge=None):
        """
        If `from_element` and `corner` are None, we compute it from this position.

        :param ep3:
        :param from_element: We compute it from this element.
        :param corner_edge: We compute it from this corner.
        :return:
        """
        if self._ee_.IS.on_periodic_boundary:
            assert from_element is not None, \
                "to compute the physical position of an edge element on periodic " \
                "boundary, we have to provide from which element you " \
                "want to compute it since it clearly will gives " \
                "different results."
            if from_element == 'any':
                # the different results do not matter; for example, when
                # we want to get value from a periodic function, the
                # location for evaluating the function also does not
                # matter.
                from_element = self._ee_.CHARACTERISTIC_element
            else:
                pass

        if from_element is None:
            i = self._ee_.CHARACTERISTIC_element
        elif from_element == 'any':
            i = self._ee_.CHARACTERISTIC_element
        else:
            i = from_element

        assert self._ee_.i in self._ee_._elements_.map[i], \
            f"edge element #{self._ee_.i} is not on mesh element {i}."

        ___ = ['WB', 'EB', 'WF', 'EF', 'NB', 'SB', 'NF', 'SF', 'NW', 'SW', 'NE', 'SE']
        if self._ee_._elements_.map[i].count(self._ee_.i) == 1:
            # this mesh element is not periodic to itself.
            corner_index = self._ee_._elements_.map[i].index(self._ee_.i)
            element_corner = ___[corner_index]
            if from_element is None: # if we do not provide `from_element` we must have this
                assert element_corner == self._ee_.CHARACTERISTIC_corner_edge

            if corner_edge is not None:
                assert element_corner == corner_edge, \
                    f"cannot compute it at provided corner edge: {corner_edge}"

        elif self._ee_._elements_.map[i].count(self._ee_.i) > 1:
            # this mesh element is periodic to itself.
            assert corner_edge is not None, f"edge element #{self._ee_.i} " \
                                            f"is on more than 1 corner-edges of element #{i} " \
                                            f"(periodic), provide corner_edge as well."
            element_corner = corner_edge
        else:
            raise Exception()

        # we will compute the physical position of this edge element from mesh element #`i` at its corner_edge `element_corner`
        assert self._ee_.i == self._ee_._elements_.map[i][___.index(element_corner)], \
            f"node element #{self._ee_.i} is not at {element_corner} of mesh element #{i}."

        ep = self._ee_._elements_.___generate_full_ep___(ep3, element_corner)
        x, y, z = self._ee_._elements_._mesh_.elements[i].coordinate_transformation.mapping(*ep)

        return x, y, z






if __name__ == '__main__':
    # mpiexec -n 6 python _3dCSCG\mesh\edge\elements\element\coordinate_transformation.py
    from objects.CSCG._3d.master import MeshGenerator
    elements = [2, 2, 2]
    # mesh = MeshGenerator('crazy_periodic', c=0.0, bounds=([0,3], [0,3], [0,3]))(elements)
    mesh = MeshGenerator('bridge_arch_cracked')(elements)
    edges = mesh.edge.elements

    for i in edges:
        edge = edges[i]