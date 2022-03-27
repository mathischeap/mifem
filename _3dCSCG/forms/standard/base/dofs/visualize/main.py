""""""


from screws.freeze.main import FrozenOnly


from _3dCSCG.forms.standard.base.dofs.visualize.matplot import _3dCSCG_SF_DOFs_VISUALIZE_matplot


class _3dCSCG_SF_DOFs_VISUALIZE(FrozenOnly):
    """"""
    def __init__(self, dofs):
        """"""
        self._dofs_ = dofs
        self._mesh_ = dofs._sf_.mesh
        self._matplot_ = None
        self._freeze_self_()

    def __call__(self, *args, **kwargs):
        return self.matplot(*args, **kwargs)

    @property
    def matplot(self):
        if self._matplot_ is None:
            self._matplot_ = _3dCSCG_SF_DOFs_VISUALIZE_matplot(self._dofs_)
        return self._matplot_
