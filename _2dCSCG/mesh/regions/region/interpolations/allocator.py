from screws.frozen import FrozenOnly
from importlib import import_module


class InterpolationSearcher(FrozenOnly):
    """ """

    def __init__(self, ID):
        """ """
        assert ID in self.___coded_interpolator___(), \
            " <InterpolationSearcher> : interpolation named '{}' is not coded.".format(ID)
        self._ID_ = ID
        cls_name = self.___coded_interpolator___()[ID]
        cls_path = self.___interpolator_path___() + '.' + ID
        self._cls_ = getattr(import_module(cls_path), cls_name)
        self._freeze_self_()

    def __call__(self, region):
        """ """
        return self._cls_(region)

    @classmethod
    def ___coded_interpolator___(cls):
        """ Update this whenever we code a new interpolator. """
        return {'transfinite': 'Transfinite',
                'crazy': 'Crazy', }

    @classmethod
    def ___interpolator_path___(cls):
        """ """
        return "_2dCSCG.mesh.regions.region.interpolations"