# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/25 6:15 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from objects.mpRfT._2d.forms.segment.base.main import mpRfT2_SegmentFormBase
from objects.mpRfT._2d.forms.segment.node.numbering.main import mpRfT2_NSgF_Numbering
from objects.mpRfT._2d.forms.segment.node.discretize.main import mpRfT2_NSgF_Discretize
from objects.mpRfT._2d.forms.segment.node.reconstruct.main import mpRfT2_NSgF_Reconstruct
from objects.mpRfT._2d.forms.segment.node.migration import mpRfT2_NSgF_Migration
from objects.mpRfT._2d.forms.segment.node.visualize import mpRfT2_NSgF_Visualize
from objects.mpRfT._2d.forms.segment.node.num import mpRfT2_NSgF_Num




class mpRfT2_NSgF(mpRfT2_SegmentFormBase):
    """"""
    def __init__(self, mesh, ndp=-1, ntype='Gauss', numbering_parameters='Naive', name='node-segment-form'):
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
        super(mpRfT2_NSgF, self).__init__(mesh, name, ndp, ntype)
        self.standard_properties.___PRIVATE_add_tag___('mpRfT2_node_segment_form')

        self._numbering_ = mpRfT2_NSgF_Numbering(self, numbering_parameters)
        self._num_ = mpRfT2_NSgF_Num(self)

        self._discretize_ = mpRfT2_NSgF_Discretize(self)
        self._reconstruct_ = mpRfT2_NSgF_Reconstruct(self)
        self._migrate_ = mpRfT2_NSgF_Migration(self)
        self._visualize_ = mpRfT2_NSgF_Visualize(self)

        self._freeze_self_()





if __name__ == "__main__":
    # mpiexec -n 4 python objects/mpRfT/_2d/forms/segment/node/main.py

    from __init__ import rfT2

    fc = rfT2.rf(100, N_range=(2,2))

    t0 = fc('nst', ndp=0, ntype='Gauss')

    mesh = t0.mesh

    import numpy as np

    from objects.mpRfT._2d.cf.scalar.main import mpRfT2_Scalar
    def p(t, x, y): return np.sin(np.pi*x) * np.sin(np.pi*y) + t
    s = mpRfT2_Scalar(mesh, p)

    t0.TW.func = s
    s.current_time = 0
    t0.discretization()
    t0.visualization()

    # print(t0.num.local_dofs)