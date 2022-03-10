
from screws.freeze.main import FrozenOnly
from importlib import import_module


class TypeWr2MetricGiver(FrozenOnly):
    """
    We use this class to distribute type object to regions.
    """
    def __init__(self, ID: str):
        assert ID in self.___coded_typeWr2Metric___(), \
            " <TypeWr2MetricGiver> : typeWr2Metric named '{}' is not coded.".format(ID)
        self._ID_ = ID
        cls_name = self.___coded_typeWr2Metric___()[ID]
        cls_path = self.___typeWr2Metric_path___() + '.' + ID
        self._itp_ = getattr(import_module(cls_path), cls_name)
        self._freeze_self_()

    def __call__(self, region):
        return self._itp_(region)

    @classmethod
    def ___coded_typeWr2Metric___(cls):
        """Update this whenever we code a new Interpolator. """
        return {'chaotic': 'Chaotic',
                'crazy': 'Crazy',
                'transfinite': 'Transfinite'}

    @classmethod
    def ___typeWr2Metric_path___(cls):
        return "_2dCSCG.mesh.domain.regions.region.types_wrt_metric"