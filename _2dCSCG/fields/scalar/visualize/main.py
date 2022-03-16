

from screws.freeze.main import FrozenOnly
from _2dCSCG.fields.scalar.visualize.matplot import _2dCSCG_ScalarField_Visualize_matplot

class _2dCSCG_ScalarField_Visualize(FrozenOnly):
    def __init__(self, cf):
        self._cf_ = cf
        self._default_ = 'matplot'
        self._mesh_ = self._cf_.mesh
        self._matplot_ = None
        self._freeze_self_()

    def __call__(self, **kwargs):
        return getattr(self, self._default_)(**kwargs)

    @property
    def matplot(self):
        """ """
        if self._matplot_ is None:
            self._matplot_ = _2dCSCG_ScalarField_Visualize_matplot(self._cf_)
        return self._matplot_
