

from screws.freeze.base import FrozenOnly


class _3dCSCG_TraceElement_VIS(FrozenOnly):
    """

    """
    def __init__(self, trace_element):
        """"""
        self._te_ = trace_element
        self._matplot_ = None
        self._default_ = 'matplot'
        self._freeze_self_()

    def __call__(self, *args, **kwargs):
        """With this method, we will only plot the trace-elements on its local characteristic
        mesh element to avoid communication between cores.

        :param kwargs:
        :return:
        """
        # we cannot use self._te_._elements_.do.illustrate_element(self._te_.i, **kwargs) because
        # that method will make use of other cores while trace elements are stored locally.
