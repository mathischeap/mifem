

from SCREWS.frozen import FrozenOnly
from _3dCSCG.form.standard.dofs.dof.visualize.matplot import _3dCSCG_SF_DOF_VISUALIZE_matplot


class _3dCSCG_SF_DOF_VISUALIZE(FrozenOnly):
    """"""
    def __init__(self, dof):
        """"""
        self._dof_ = dof
        self._mesh_ = dof._sf_.mesh
        self._matplot_ = _3dCSCG_SF_DOF_VISUALIZE_matplot(dof)
        self._freeze_self_()

    def __call__(self, *args, **kwargs):
        return self.matplot(*args, **kwargs)

    @property
    def matplot(self):
        return self._matplot_
