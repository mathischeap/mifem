



from screws.freeze.main import FrozenOnly
from typing import Dict
from objects.CSCG.base.ADF.base.BC.partial_cochain.partial_dofs.include_from import \
    _ADF_PartialDofs_Include_from_
from objects.CSCG.base.ADF.base.BC.partial_cochain.partial_dofs.interpretation.main import \
    _ADF_PartialDofs_Interpretation_



class PartialDofs(FrozenOnly):
    """"""
    def __init__(self, adf):
        assert 'CSCG_form' in adf.prime.standard_properties.tags
        self._adf_ = adf
        self._mesh_ = adf.mesh
        self._dofs_: Dict[int, list] = dict()
        self._include_ = _ADF_PartialDofs_Include_from_(self)
        self._interpretation_ = _ADF_PartialDofs_Interpretation_(self)
        self._freeze_self_()

    @property
    def dofs(self):
        """Return a dict. The keys are mesh element numbers. And values
        are lists of indicators indicating which local dofs.

        indicators:
            type-1: '1-' + 'N', 'S', 'W', 'E', 'B', 'F' for 3D or 'U', 'D', 'L', 'R'
                for 2D CSCG mesh --- the dofs on mesh element side.

        """
        return self._dofs_

    @property
    def include(self):
        """Methods used to include local dofs."""
        return self._include_

    @property
    def interpreted_as(self):
        return self._interpretation_
