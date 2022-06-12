# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/06/01 5:02 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from objects.mpRfT._2d.forms.segment.base.main import mpRfT2_SegmentFormBase

from objects.mpRfT._2d.forms.segment.edge.numbering.main import mpRfT2_ESgF_Numbering
from objects.mpRfT._2d.forms.segment.edge.discretize.main import mpRfT2_ESgF_Discretize
from objects.mpRfT._2d.forms.segment.edge.reconstruct.main import mpRfT2_ESgF_Reconstruct
from objects.mpRfT._2d.forms.segment.edge.visualize import mpRfT2_ESgF_Visualize
from objects.mpRfT._2d.forms.segment.edge.num import mpRfT2_ESgF_Num




class mpRfT2_ESgF(mpRfT2_SegmentFormBase):
    """"""
    def __init__(self, mesh, ndp=0, ntype='Lobatto', numbering_parameters='Naive', name='edge-segment-form'):
        """

        Parameters
        ----------
        mesh
        ndp : int
            node degree plus.
        ntype : str
            node type.
        numbering_parameters :
        name :

        """
        super(mpRfT2_ESgF, self).__init__(mesh, name, ndp, ntype)
        self.standard_properties.___PRIVATE_add_tag___('mpRfT2_edge_segment_form')

        self._numbering_ = mpRfT2_ESgF_Numbering(self, numbering_parameters)
        self._num_ = mpRfT2_ESgF_Num(self)

        self._discretize_ = mpRfT2_ESgF_Discretize(self)
        self._reconstruct_ = mpRfT2_ESgF_Reconstruct(self)
        self._visualize_ = mpRfT2_ESgF_Visualize(self)

        self._freeze_self_()



if __name__ == "__main__":
    # mpiexec -n 4 python objects/mpRfT/_2d/forms/segment/edge/main.py

    from __init__ import rfT2

    fc = rfT2.rf(100, N_range=(2,2))

    t = fc('est', ndp=0, ntype='Lobatto')

    mesh = t.mesh

    # import numpy as np
    #
    # from objects.mpRfT._2d.cf.scalar.main import mpRfT2_Scalar
    # def p(t, x, y): return np.sin(np.pi*x) * np.sin(np.pi*y) + t
    # s = mpRfT2_Scalar(mesh, p)
    #
    # t0.TW.func = s
    # s.current_time = 0
    # t0.discretize()
    # t0.visualize()
