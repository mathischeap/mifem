# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
from screws.freeze.main import FrozenOnly


from _2dCSCG.forms.standard.base.visualize.matplot import _2dCSCG_FormVisualize_Matplot

class _2dCSCG_FormVisualize(FrozenOnly):
    """The visualization property/component of standard forms."""
    def __init__(self, sf):
        self._sf_ = sf
        self._defaultPlot_ = 'matplot'
        self._mesh_ = sf.mesh
        self._matplot_ = _2dCSCG_FormVisualize_Matplot(sf)
        self._freeze_self_()

    def __call__(self, **kwargs):
        """When this object is called, we call the default visualizing method: ``tecplot``."""
        return getattr(self, self._defaultPlot_)(**kwargs)

    @property
    def matplot(self):
        return self._matplot_

