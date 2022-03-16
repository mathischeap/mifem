


from screws.freeze.main import FrozenOnly
from _2dCSCG.APP.exact_solution.visualize.matplot import _2dCSCG_ES_VIS_Matplot

class ExactSolution_Visualize(FrozenOnly):
    """Here we try to visualize the ES from a higher level than individual variables, for those,
    we can call the visualizing method of themselves.
    """
    def __init__(self, es):
        self._es_ = es
        self._matplot_ = _2dCSCG_ES_VIS_Matplot(es)
        self._freeze_self_()


    def matplot(self):
        """"""
        return self._matplot_