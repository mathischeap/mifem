# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/08 10:29 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
import numpy as np

class SegmentCT(FrozenOnly):
    """"""

    def __init__(self, sg):
        """"""
        self._sg_ = sg
        self._freeze_self_()


    def mapping(self, xi):
        """from [-1,1] to this local segment."""

        sg = self._sg_
        name = sg.where.__class__.__name__
        if name == '_2nCSCG_Lv0TraceElement':
            return self._Pr_mapping_2nCSCG_Lv0TraceElement(xi)
        elif name == '_2nCSCG_RF2_MeshCell':
            return self._Pr_mapping_2nCSCG_RF2_MeshCell(xi)
        else:
            raise NotImplementedError(
                f"segment mapping not implemented for sg.where={sg.where.__class__.__name__}")


    def _Pr_mapping_2nCSCG_Lv0TraceElement(self, xi):
        """"""
        sg = self._sg_
        o, d = sg.origin_and_delta
        xi = o + (xi + 1) * d / 2
        return sg.where._te_.coordinate_transformation.mapping(xi)


    def _Pr_mapping_2nCSCG_RF2_MeshCell(self, xi):
        """"""
        sg = self._sg_
        o, d = sg.origin_and_delta
        ox, oy = o

        if isinstance(sg.xSignature, str):
            xi = ox + (xi + 1) * d / 2
            et = oy * np.ones_like(xi)
        else:
            et = oy + (xi + 1) * d / 2
            xi = ox * np.ones_like(et)

        return sg.where.coordinate_transformation.mapping(xi, et)



if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
