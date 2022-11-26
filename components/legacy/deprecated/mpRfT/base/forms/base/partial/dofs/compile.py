# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/06/13 9:11 PM
"""
import sys

if './' not in sys.path: sys.path.append('/')

from components.freeze.base import FrozenOnly


class mpRfT_Form_PartialDofs_Compile(FrozenOnly):
    """"""

    def __init__(self, pd):
        """"""
        self._f_ = pd._f_
        self._mesh_ = pd._mesh_
        self._pd_ = pd
        self._freeze_self_()

    @property
    def local(self):
        """The partial dofs are compiled as root-cell-wise local-dofs.

        Return a dict whose keys are the repr of the root-cells, and values are the local dofs.
        """
        mpRfT2_edges = ['U', 'D', 'L', 'R']
        indicators = self._pd_._indicators_

        DOFs =  dict()
        for rc_rp in indicators:
            DOFs[rc_rp] = list()
            indS = indicators[rc_rp]
            for indicator in indS:

                #------------------------------------------------------------------------------
                if indicator in mpRfT2_edges:
                    dofs = self._f_.numbering.do.find.local_numbering_of_dofs_on(rc_rp, indicator)

                #------------------------------------------------------------------------------
                else:
                    raise NotImplementedError(f"Not implemented for indicator={indicator}")
                #==============================================================================

                DOFs[rc_rp].extend(dofs)

        return DOFs







if __name__ == "__main__":
    # mpiexec -n 4 python objects/mpRfT/base/forms/base/partial/dofs/compile.py
    from __init__ import rfT2

    fc = rfT2.rf(100, N_range=(2,3), mesh_pool='crazy')

    f1 = fc('1-f-o')
    f1.BC.valid_boundaries = ['Upper', 'Left']
    fpd = f1.BC.partial_dofs

    t1 = fc('nst')
    t1.BC.valid_boundaries = ['Upper', 'Left']
    tpd = t1.BC.partial_dofs

    dofs_f = fpd.compile.local
    dofs_t = tpd.compile.local

    print(dofs_f)
    print(dofs_t)