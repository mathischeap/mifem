


from screws.freeze.main import FrozenOnly
from objects.CSCG._3d.forms.standard.base.visualize.matplot import _3dCSCG_standard_form_Matplot
from objects.CSCG._3d.forms.standard.base.visualize.tecplot_ import _3dCSCG_standard_form_Tecplot




class _3dCSCG_FormVisualize(FrozenOnly):
    """The visualization property/component of standard forms."""
    def __init__(self, sf):
        assert '3dCSCG_standard_form' in sf.standard_properties.tags
        self._sf_ = sf
        self._defaultPlot_ = 'tecplot'
        self._tecplot_= _3dCSCG_standard_form_Tecplot(sf)
        self._matplot_= _3dCSCG_standard_form_Matplot(sf)
        self._freeze_self_()

    def __call__(self, **kwargs):
        """When this object is called, we call the default visualizing method: ``tecplot``."""
        return getattr(self, self._defaultPlot_)(**kwargs)

    @property
    def tecplot(self):
        return self._tecplot_

    @property
    def matplot(self):
        return self._matplot_