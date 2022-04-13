


from screws.freeze.main import FrozenOnly
from objects.CSCG._3d.forms.standard.base.visualize.tecplot_ import _3dCSCG_standard_form_Tecplot




class _3dCSCG_FormVisualize(FrozenOnly):
    """The visualization property/component of standard forms."""
    def __init__(self, sf):
        assert '3dCSCG_standard_form' in sf.standard_properties.tags
        self._tecplot_= _3dCSCG_standard_form_Tecplot(sf)

    @property
    def tecplot(self):
        return self._tecplot_
