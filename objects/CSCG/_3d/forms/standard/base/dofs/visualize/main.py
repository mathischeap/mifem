# -*- coding: utf-8 -*-
""""""
from screws.freeze.main import FrozenOnly

from objects.CSCG._3d.forms.standard.base.dofs.visualize.matplot._0sf import _3dCSCG_S0F_DOFs_Matplot
from objects.CSCG._3d.forms.standard.base.dofs.visualize.matplot._1sf import _3dCSCG_S1F_DOFs_Matplot
from objects.CSCG._3d.forms.standard.base.dofs.visualize.matplot._2sf import _3dCSCG_S2F_DOFs_Matplot
from objects.CSCG._3d.forms.standard.base.dofs.visualize.matplot._3sf import _3dCSCG_S3F_DOFs_Matplot


class _3dCSCG_SF_DOFs_VISUALIZE(FrozenOnly):
    """"""

    def __init__(self, dofs):
        """"""
        self._dofs_ = dofs
        self._matplot_ = None
        self._freeze_self_()

    def __call__(self, *args, **kwargs):
        return self.matplot(*args, **kwargs)

    @property
    def matplot(self):
        if self._matplot_ is None:

            if self._dofs_._sf_.k == 0:
                self._matplot_ = _3dCSCG_S0F_DOFs_Matplot(self._dofs_)

            elif self._dofs_._sf_.k == 1:
                self._matplot_ = _3dCSCG_S1F_DOFs_Matplot(self._dofs_)

            elif self._dofs_._sf_.k == 2:
                self._matplot_ = _3dCSCG_S2F_DOFs_Matplot(self._dofs_)

            elif self._dofs_._sf_.k == 3:
                self._matplot_ = _3dCSCG_S3F_DOFs_Matplot(self._dofs_)

            else:
                raise Exception()

        return self._matplot_