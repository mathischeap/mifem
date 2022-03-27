# -*- coding: utf-8 -*-
"""
The trace (face) elements of a mesh.
"""

import sys
if './' not in sys.path: sys.path.append('../')

from root.config.main import *
from screws.freeze.main import FrozenOnly

from _3dCSCG.mesh.trace.elements.main import _3dCSCG_Trace_Elements


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

        AvQ = cOmm.gather(AvQ, root=mAster_rank)
        WorstQ = cOmm.gather(WorstQ, root=mAster_rank)
        BestQ = cOmm.gather(BestQ, root=mAster_rank)

        if rAnk == mAster_rank:
            AvQ = [i for i in AvQ if i]
            WorstQ = [i for i in WorstQ if i]
            BestQ = [i for i in BestQ if i]
            AvQ = sum(AvQ) / len(AvQ)
            WorstQ = min(WorstQ)
            BestQ = max(BestQ)
        else:
            pass
        AvQ = cOmm.bcast(AvQ, root=mAster_rank)
        WorstQ = cOmm.bcast(WorstQ, root=mAster_rank)
        BestQ = cOmm.bcast(BestQ, root=mAster_rank)

        return {'average quality': AvQ,
                'worst quality': WorstQ,
                'best quality': BestQ}






if __name__ == '__main__':
    # mpiexec -n 12 python _3dCSCG\mesh\trace.py
    from _3dCSCG.master import MeshGenerator
    elements = [3, 4, 2]
    mesh = MeshGenerator('crazy_periodic', c=0.3, bounds=([0,1], [0,1], [0,1]))(elements)
    mesh.trace.elements.SELFCHECK.outward_unit_normal_vector()
    Q = mesh.trace.elements.quality
    print(mesh.quality)
    print(mesh.trace.quality)