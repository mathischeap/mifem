# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
from screws.freeze.main import FrozenOnly


class _2dCSCG_TraceVisualize(FrozenOnly):
    """The visualization property/component of standard forms."""
    def __init__(self, tf):
        self._tf_ = tf
        self._defaultPlot_ = 'matplot'
        self._freeze_self_()

    def __call__(self, **kwargs):
        """When this object is called, we call the default visualizing method: ``tecplot``."""
        return getattr(self, self._defaultPlot_)(**kwargs)

    def matplot(self, **kwargs):
        """Visualize the standard form using ``tecplot``."""
        getattr(self, f"_matplot_{self._tf_.__class__.__name__}")(**kwargs)


    def _matplot_1Trace_Outer(self):
        """Tecplot for 2-trace-form."""
        print('_matplot_1Trace_Outer')