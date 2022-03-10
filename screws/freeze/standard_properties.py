
import sentry_sdk
from root.config.main import seNtry_on
if seNtry_on: sentry_sdk.init("https://79fb951c3ea7457c8c3f47c3dfb026ce@sentry.io/1458280")

from screws.exceptions import StatisticError
from screws.exceptions import ParametersError
from screws.freeze.inheriting.frozen_only import FrozenOnly


class StandardProperties(FrozenOnly):
    def __init__(self, obj):
        self.___mark___ = None
        self.___name___ = None
        self.___stamp___ = None
        self.___tags___ = tuple()
        self.___obj___ = obj
        self._freeze_self_()

    @property
    def name(self):
        """
        ``name``, ``mark``, ``stamp`` are just three standard reserved property.

        :return:
        """
        return self.___name___

    @name.setter
    def name(self, name):
        if name is None:
            self.___name___ = None
        else:
            assert isinstance(name, str), " name needs to be a string."
            self.___name___ = name


    @property
    def mark(self):
        """
        ``name``, ``mark``, ``stamp`` are just three standard reserved property.

        :return:
        """
        return self.___mark___

    @mark.setter
    def mark(self, mark):
        self.___mark___ = mark


    @property
    def stamp(self):
        """
        ``name``, ``mark``, ``stamp`` are just three standard reserved property.

        :return:
        """
        return self.___stamp___

    @stamp.setter
    def stamp(self, stamp):
        self.___stamp___ = stamp


    @property
    def parameters(self):
        """
        Parameters must be enough for uniquely define an object.

        So, we can recover an object singly from its parameters.
        """

        try:
            if not hasattr(self.___obj___, '___parameters___'):
                raise ParametersError
        except ParametersError as PE:
            if seNtry_on:
                sentry_sdk.capture_exception(PE)
            raise ParametersError(
                " <Parameters> : %r has no ___parameters___ property" % self.___obj___)
        return self.___obj___.___parameters___

    @property
    def statistic(self):
        """
        statistic of an object returns a dict contain some statistic info of the object.

        It can not be used to recover an object, neither does not represent an object uniquely.

        :return:
        """
        try:
            if not hasattr(self.___obj___, '___statistic___'):
                raise StatisticError
        except StatisticError as SE:
            if seNtry_on:
                sentry_sdk.capture_exception(SE)
            raise StatisticError(
                " <Statistic> : %r has no ___statistic___ property" % self.___obj___)
        return self.___obj___.___statistic___

    @property
    def tags(self):
        """
        tags are just tags of an object, we use them to just tag it.

        A tag can not be removed.
        """
        return self.___tags___

    def ___PRIVATE_add_tag___(self, tag):
        """"""
        self.___tags___ += (tag,)