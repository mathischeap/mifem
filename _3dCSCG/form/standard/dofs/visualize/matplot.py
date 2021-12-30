
from SCREWS.frozen import FrozenOnly



class _3dCSCG_SF_DOFs_VISUALIZE_matplot(FrozenOnly):
    """"""
    def __init__(self, dofs):
        """"""
        self._dofs_ = dofs
        self._sf_ = dofs._sf_
        self._mesh_ = dofs._sf_.mesh
        self._freeze_self_()


    def __call__(self, *args, **kwargs):
        """Default visualizer"""
        return getattr(self, f"___PRIVATE_matplot_local_dofs_{self._sf_.k}form___")(*args, **kwargs)


    def ___PRIVATE_matplot_local_dofs_0form___(self):
        """We will call all cores to do plots for themselves all together."""
        raise NotImplementedError()

    def ___PRIVATE_matplot_local_dofs_1form___(self):
        """We will call all cores to do plots for themselves all together."""
        raise NotImplementedError()

    def ___PRIVATE_matplot_local_dofs_2form___(self):
        """We will call all cores to do plots for themselves all together."""
        raise NotImplementedError()

    def ___PRIVATE_matplot_local_dofs_3form___(self):
        """We will call all cores to do plots for themselves all together."""
        raise NotImplementedError()



