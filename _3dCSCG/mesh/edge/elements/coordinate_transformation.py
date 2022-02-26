

import sys
if './' not in sys.path: sys.path.append('./')

from screws.frozen import FrozenOnly


class _3dCSCG_Edge_Elements_CT(FrozenOnly):
    """"""
    def __init__(self, elements):
        """"""
        self._elements_ = elements
        self._freeze_self_()






if __name__ == '__main__':
    # mpiexec -n 12 python _3dCSCG\mesh\edge\elements\coordinate_transformation.py
    from _3dCSCG.main import MeshGenerator
    elements = [2, 2, 2]
    # mesh = MeshGenerator('crazy_periodic', c=0.0, bounds=([0,3], [0,3], [0,3]))(elements)
    mesh = MeshGenerator('bridge_arch_cracked')(elements)
    edges = mesh.edge.elements
