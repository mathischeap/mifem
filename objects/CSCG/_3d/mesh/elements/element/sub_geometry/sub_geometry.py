
"""
Sub-geometries of an element. Like we can pick a point, a slice or a volume from an element.

"""

from screws.freeze.main import FrozenOnly
from objects.CSCG._3d.mesh.elements.element.sub_geometry.helpers.perpendicular_slice import \
    ElementPerpendicularSlice


class ElementSubGeometry(FrozenOnly):
    """"""

    def __init__(self, element):
        """"""
        self._element_ = element
        self._freeze_self_()


    def make_a_perpendicular_slice_object_on(self, xi=None, eta=None, sigma=None):
        """Only one of ``xi``, ``eta`` and ``sigma`` can be a float in :math:`[-1,1]`, the other two must be None.

        :param xi:
        :param eta:
        :param sigma:
        :return:
        """
        return ElementPerpendicularSlice(self._element_, xi=xi, eta=eta, sigma=sigma)