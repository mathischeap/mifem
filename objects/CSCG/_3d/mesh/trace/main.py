# -*- coding: utf-8 -*-
"""
The trace (face) elements of a mesh.
"""

import sys
if './' not in sys.path: sys.path.append('./')

from root.config.main import *
from components.freeze.main import FrozenOnly

from objects.CSCG._3d.mesh.trace.elements.main import _3dCSCG_Trace_Elements


class _3dCSCG_Trace(FrozenOnly):
    def __init__(self, mesh):
        self._mesh_ = mesh
        self._elements_ = _3dCSCG_Trace_Elements(self) # please initialize it here!
        self.___PRIVATE_reset_cache___()
        self._freeze_self_()

    def ___PRIVATE_reset_cache___(self):
        pass

    @property
    def elements(self):
        return self._elements_

    @property
    def quality(self):
        _ = self.elements.quality
        AvQ, WorstQ, BestQ = self.elements._AvQ_, \
                             self.elements._WorstQ_, \
                             self.elements._BestQ_

        AvQ = COMM.gather(AvQ, root=MASTER_RANK)
        WorstQ = COMM.gather(WorstQ, root=MASTER_RANK)
        BestQ = COMM.gather(BestQ, root=MASTER_RANK)

        if RANK == MASTER_RANK:
            AvQ = [i for i in AvQ if i]
            WorstQ = [i for i in WorstQ if i]
            BestQ = [i for i in BestQ if i]
            AvQ = sum(AvQ) / len(AvQ)
            WorstQ = min(WorstQ)
            BestQ = max(BestQ)
        else:
            pass
        AvQ = COMM.bcast(AvQ, root=MASTER_RANK)
        WorstQ = COMM.bcast(WorstQ, root=MASTER_RANK)
        BestQ = COMM.bcast(BestQ, root=MASTER_RANK)

        return {'average quality': AvQ,
                'worst quality': WorstQ,
                'best quality': BestQ}






if __name__ == '__main__':
    # mpiexec -n 6 python objects/CSCG/_3d/mesh/trace/main.py
    from objects.CSCG._3d.master import MeshGenerator
    elements = [3, 4, 2]
    mesh = MeshGenerator('crazy_periodic', c=0.3, bounds=([0,1], [0,1], [0,1]))(elements)
    mesh.trace.elements.selfcheck.outward_unit_normal_vector()
    Q = mesh.trace.elements.quality
    print(mesh.quality)
    print(mesh.trace.quality)