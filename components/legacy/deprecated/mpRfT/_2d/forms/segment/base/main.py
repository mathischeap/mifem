# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/25 6:12 PM
"""
import sys

if './' not in sys.path: sys.path.append('/')

from objects.mpRfT._2d.forms.base import mpRfT2_FormBase

from objects.mpRfT._2d.forms.segment.base.cochain import mpRfT2_SgF_Cochain
from objects.mpRfT._2d.forms.segment.base.IS import mpRfT2_SgF_IS
from objects.mpRfT._2d.forms.segment.base.N import mpRfT2_SgF_N


class mpRfT2_SegmentFormBase(mpRfT2_FormBase):
    """"""

    def __init__(self, mesh, name, ndp, ntype):
        """

        Parameters
        ----------
        mesh
        name
        ndp
        ntype
        """
        super(mpRfT2_SegmentFormBase, self).__init__(mesh, name)
        self.standard_properties.___PRIVATE_add_tag___('mpRfT2_segment_form')

        self.ndp = ndp
        self._ntype_= ntype
        self._cochain_ = mpRfT2_SgF_Cochain(self)
        self._IS_ = mpRfT2_SgF_IS(self)
        self._N_ = mpRfT2_SgF_N(self)

        self._numbering_ = None
        self._num_ = None

        self._discretize_ = None
        self._reconstruct_ = None
        self._migrate_ = None
        self._visualize_ = None

    #-------- must have methods --------------------------------------------------------------
    def ___Pr_check_analytic_expression___(self, func):
        """"""
        assert func.mesh is self.mesh

        if func.__class__.__name__ == 'mpRfT2_Scalar':
            assert func.ftype in ('standard',), \
                f"mpRfT2-trace-form FUNC do not accept func mpRfT2_Scalar of ftype {func.ftype}."
        else:
            raise Exception(f"mpRfT2-trace-form FUNC do not accept func {func.__class__}")


    def ___Pr_check_BC_analytic_expression___(self, ae):
        assert ae.mesh is self.mesh

        if ae.__class__.__name__ == 'mpRfT2_Scalar':
            assert ae.ftype in ('standard',), \
                f"mpRfT2-trace-form BC do not accept func mpRfT2_Scalar of ftype {ae.ftype}."
        else:
            raise Exception(f"mpRfT2-trace-form BC do not accept func {ae.__class__}")


    #-------------------------------------------------------------------------------------------
    @property
    def ndp(self):
        return self._ndp_

    @ndp.setter
    def ndp(self, ndp):
        """

        Parameters
        ----------
        ndp :
            N degree plus.

        Returns
        -------

        """
        if isinstance(ndp, int):
            self._ndp_ = ndp
        else:
            raise NotImplementedError(f"Cannot handle ndp = {ndp}.")

    @property
    def N(self):
        """"""
        return self._N_

    @property
    def ntype(self):
        return self._ntype_

    @property
    def cochain(self):
        return self._cochain_

    @property
    def IS(self):
        return self._IS_

    #-----------------------------------------------------------------------------------------
    @property
    def numbering(self):
        return self._numbering_

    @property
    def num(self):
        return self._num_

    #-----------------------------------------------------------------------------------------
    @property
    def discretization(self):
        return self._discretize_

    @property
    def reconstruction(self):
        return self._reconstruct_

    @property
    def migration(self):
        return self._migrate_

    @property
    def visualization(self):
        return self._visualize_






if __name__ == "__main__":
    # mpiexec -n 4 python objects/mpRfT/_2d/forms/segment/base/main.py
    pass
