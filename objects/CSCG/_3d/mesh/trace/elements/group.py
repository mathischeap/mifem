

import sys
if './' not in sys.path: sys.path.append('/')

from screws.freeze.main import FrozenOnly



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




if __name__ == '__main__':
    # mpiexec -n 12 python _3dCSCG\mesh\trace\elements\group.py
    from objects.CSCG._3d.master import MeshGenerator
    elements = [3, 4, 2]
    mesh = MeshGenerator('crazy_periodic', c=0.3, bounds=([0,1], [0,1], [0,1]))(elements)
    mesh.trace.elements.SELFCHECK.outward_unit_normal_vector()
    Q = mesh.trace.elements.quality

    mesh.trace.elements.do.illustrate_trace_element(1)