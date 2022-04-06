from screws.freeze.main import FrozenOnly


from objects.CSCG.base.ADF.base.BC.partial_cochain.include_from import _ADF_PartialCochain_Include_from_
from objects.CSCG.base.ADF.base.BC.partial_cochain.interpretation import _ADF_PartialCochain_Interpretation_
from objects.CSCG.base.ADF.base.BC.partial_cochain.partial_dofs.main import PartialDofs

class PartialCochain(FrozenOnly):
    """"""
    def __init__(self, adf):
        assert 'CSCG_form' in adf.prime.standard_properties.tags
        self._adf_ = adf
        self._mesh_ = adf.mesh
        self._dofs_ = PartialDofs(adf)
        self._cochain_ = dict()
        self._include_ = _ADF_PartialCochain_Include_from_(self)
        self._interpretation_ = _ADF_PartialCochain_Interpretation_(self)
        self._freeze_self_()

    @property
    def dofs(self):
        """The dofs are included in a PartialDofs instance."""
        return self._dofs_

    @property
    def cochain(self):
        """A dict: keys are mesh element numbers, values are the local cochains. So it is corresponding to
        `self.dofs.interpreted_as.local_dofs`.
        """
        return self._cochain_

    @property
    def include(self):
        return self._include_

    @property
    def interpreted_as(self):
        return self._interpretation_