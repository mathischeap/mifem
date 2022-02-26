


from screws.frozen import FrozenOnly


from _3dCSCG.spaces.base.visualize.matplot import _3dCSC_Space_Visualize_Matplot



class _3dCSC_Space_Visualize(FrozenOnly):
    """"""
    def __init__(self, space):
        self._space_ = space
        self._matplot_ = _3dCSC_Space_Visualize_Matplot(space)
        self._freeze_self_()

    def __call__(self, *args, **kwargs):
        """Call the default visualizer."""
        return self.matplot(*args, **kwargs)

    @property
    def matplot(self):
        return self._matplot_