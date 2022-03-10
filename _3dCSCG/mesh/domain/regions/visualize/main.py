"""

"""

from screws.freeze.main import FrozenOnly
from _3dCSCG.mesh.domain.regions.visualize.matplot import _3dCSCG_Regions_Visualize_Matplot_


class _3dCSCG_Regions_Visualize(FrozenOnly):
    def __init__(self, regions):
        """ """
        assert regions.__class__.__name__ == 'Regions'
        self._regions_ = regions
        self._matplot_ = _3dCSCG_Regions_Visualize_Matplot_(self)
        self._freeze_self_()

    def __call__(self, *args, **kwargs):
        return self.matplot(*args, **kwargs)

    @property
    def matplot(self):
        """"""
        return self._matplot_