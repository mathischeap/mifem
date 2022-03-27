"""
We can export a form to multiple formats and data structures.


"""
from screws.freeze.main import FrozenOnly
from objects.CSCG._3d.forms.standard.base.export.field import _3dCSC_SF_Export_Field

class _3dCSC_Standard_Form_Export(FrozenOnly):
    """"""

    def __init__(self, sf):
        assert '3dCSCG_standard_form' in sf.standard_properties.tags
        self._sf_ = sf
        self._field_ = _3dCSC_SF_Export_Field(sf)

        self._freeze_self_()


    @property
    def field(self):
        return self._field_