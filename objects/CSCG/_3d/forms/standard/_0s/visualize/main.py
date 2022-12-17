# -*- coding: utf-8 -*-
from objects.CSCG._3d.forms.standard._0s.visualize.matplot import _3dCSCG_S0F_VISUALIZE_Matplot
from objects.CSCG._3d.forms.standard.base.visualize.main import _3dCSCG_FormVisualize


class _3dCSCG_S0F_VISUALIZE(_3dCSCG_FormVisualize):
    """"""
    def __init__(self, sf):
        """"""
        super(_3dCSCG_S0F_VISUALIZE, self).__init__(sf)
        self._sf_ = sf
        self._matplot_ = None
        self._freeze_self_()

    def __call__(self, *args, **kwargs):
        return self.matplot(*args, **kwargs)

    @property
    def matplot(self):
        if self._matplot_ is None:
            self._matplot_ = _3dCSCG_S0F_VISUALIZE_Matplot(self._sf_)
        return self._matplot_
