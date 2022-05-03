
from screws.exceptions import FrozenError
from abc import ABC


class FrozenOnly(ABC):

    def __setattr__(self, key, value):
        """"""
        if self.___isfrozen___ and key not in dir(self):
            raise FrozenError(" <FrozenClass> : %r is frozen. CANNOT define new attributes." % self)
        object.__setattr__(self, key, value)

    def _freeze_self_(self):
        """Freeze self, can add no more new attributes. """
        if hasattr(self, '___PreFrozenChecker___'):
            self.___PreFrozenChecker___()
        else:
            pass
        self.___ISFROZEN___ = True

    def _melt_self_(self):
        """Melt self, so  we can define new attributes."""
        self.___ISFROZEN___ = False

    @property
    def _frozen_(self):
        """Return the status of the form, frozen (True) or melt (False)?"""
        return self.___isfrozen___

    @property
    def ___isfrozen___(self):
        """"""
        try:
            return self.___ISFROZEN___
        except AttributeError:
            object.__setattr__(self, '___ISFROZEN___', False)
            return self.___ISFROZEN___