# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/21 4:34 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from tools.linear_algebra.elementwise_cache.objects.sparse_matrix.main import EWC_SparseMatrix
from objects.miUsGrid.triangular.forms.standard.base.do import miUs_Triangular_SF_Do
from objects.miUsGrid.triangular.forms.standard._0.base.do.helpers._0x1_ip_1 import ___0_x_1__ip__1___



class miUs_Triangular_S0F_Do(miUs_Triangular_SF_Do):
    """"""

    def __init__(self, sf):
        """"""
        super(miUs_Triangular_S0F_Do, self).__init__(sf)
        self._freeze_self_()


    def cross_inner(self, u, e, output='MDM'):
        """Let w denote self. Here we do (w x u, e).

        Parameters
        ----------
        u
        e
        output

        Returns
        -------

        """
        assert u.IS.standard_form and e.IS.standard_form, f"u and e must be standard form."
        assert u.mesh == self._sf_.mesh, f"u.mesh does not match."
        assert e.mesh == self._sf_.mesh, f"e.mesh does not match."

        if u.k == e.k == 1:
            DG = ___0_x_1__ip__1___(self._sf_, u, e,)
            if output == 'MDM':
                return DG.MDM
            elif output == '2-M-1':
                return EWC_SparseMatrix(self._sf_.mesh, DG.___2_M_1___, 'no_cache')
            elif output == '2-M-0':
                return EWC_SparseMatrix(self._sf_.mesh, DG.___2_M_0___, 'no_cache')
            else:
                raise NotImplementedError()
        else:
            raise NotImplementedError()




if __name__ == "__main__":
    # mpiexec -n 4 python objects/miUsGrid/triangular/forms/standard/_0/base/do/main.py
    from __init__ import miTri
    fc = miTri.form('st8', 2)
    w = fc('0-f-o')
    u = fc('1-f-o')
    e = fc('1-f-o')

    mdm = w.do.cross_inner(u, e)
    for e in mdm:
        print(e, mdm[e].shape)
