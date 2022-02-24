


import sys
if './' not in sys.path: sys.path.append('./')

from SCREWS.frozen import FrozenOnly


class _3dCSCG_Trace_Elements_DO(FrozenOnly):
    """We find some specific groups of elements."""

    def __init__(self, trace_elements):
        self._elements_ = trace_elements
        self._freeze_self_()

    def illustrate_trace_element(self, *args, **kwargs):
        return self._elements_.___DO_illustrate_trace_element___(*args, **kwargs)






if __name__ == '__main__':
    # mpiexec -n 12 python _3dCSCG\mesh\trace\elements\DO.py
    from _3dCSCG.main import MeshGenerator
    elements = [3, 4, 2]
    mesh = MeshGenerator('crazy_periodic', c=0.3, bounds=([0,1], [0,1], [0,1]))(elements)
    mesh.trace.elements.SELFCHECK.outward_unit_normal_vector()
    Q = mesh.trace.elements.quality

    mesh.trace.elements.DO.illustrate_trace_element(1)