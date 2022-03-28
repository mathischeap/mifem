



from screws.freeze.main import FrozenOnly
from objects.CSCG.base.forms.base.BC.partial_cochain.partial_dofs.interpretation.local import _PartialDofs_Interpretation_Local_
from objects.CSCG.base.forms.base.BC.partial_cochain.partial_dofs.interpretation.GLOBAL import _PartialDofs_Interpretation_Globe_




class _PartialDofs_Interpretation_(FrozenOnly):
    """A wrapper of local dof interpretation methods."""
    def __init__(self, pd):
        self._pd_ = pd
        self._mesh_ = pd._mesh_
        self._form_ = pd._form_
        self._local_dofs_ = None
        self._globe_dofs_ = None
        self._freeze_self_()

    @property
    def local_dofs(self):
        if self._local_dofs_ is None:
            self._local_dofs_ = _PartialDofs_Interpretation_Local_(self)
        return self._local_dofs_

    @property
    def globe_dofs(self):
        if self._globe_dofs_ is None:
            self._globe_dofs_ = _PartialDofs_Interpretation_Globe_(self)
        return self._globe_dofs_


