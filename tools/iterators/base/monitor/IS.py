# -*- coding: utf-8 -*-


from components.freeze.main import FrozenOnly

class IteratorMonitorIS(FrozenOnly):
    """"""
    def __init__(self, monitor):
        """"""
        self._monitor_ = monitor
        self._freeze_self_()

    @property
    def open(self):
        """(bool) Return ``True`` if the iterator is open (no max_steps, stop when shut_down)."""
        return self._monitor_._isOpen_