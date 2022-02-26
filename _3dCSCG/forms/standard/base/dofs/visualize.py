""""""


from screws.frozen import FrozenOnly





class _3dCSCG_SF_DOFs_VISUALIZE(FrozenOnly):
    """"""
    def __init__(self, dofs):
        """"""
        self._dofs_ = dofs
        self._mesh_ = dofs._sf_.mesh
        self._matplot_ = _3dCSCG_SF_DOFs_VISUALIZE_matplot(dofs)
        self._freeze_self_()

    def __call__(self, *args, **kwargs):
        return self.matplot(*args, **kwargs)

    @property
    def matplot(self):
        return self._matplot_



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



