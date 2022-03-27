from screws.freeze.base import FrozenOnly


class _2dCSCG_Outer0Trace_discretize(FrozenOnly):
    """"""
    def __init__(self, tf):
        self._tf_ = tf
        self._freeze_self_()

    def __call__(self, update_cochain=True, **kwargs):
        """
        Do the discretization.

        :param bool update_cochain: Whether we update the cochain if the trace form.
        :param kwargs: Keywords arguments to be passed to particular discretization schemes.
        :return: The cochain corresponding to the particular discretization scheme.
        """
        if self._tf_.func.ftype == 'standard':
            raise NotImplementedError()
        else:
            raise NotImplementedError()