# -*- coding: utf-8 -*-
from components.exceptions import StatisticError
from components.exceptions import ParametersError
from components.freeze.base import FrozenOnly


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
        ``name``, ``mark``, ``stamp`` ``signature`` are just standard reserved properties.

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
    def signature(self):
        """
        ``name``, ``mark``, ``stamp`` ``signature`` are just standard reserved properties.

        :return:
        """
        return self.___signature___

    @signature.setter
    def signature(self, signature):
        if signature is None:
            self.___signature___ = None
        else:
            assert isinstance(signature, str), " signature needs to be a string."
            # noinspection PyAttributeOutsideInit
            self.___signature___ = signature

    @property
    def mark(self):
        """
        ``name``, ``mark``, ``stamp`` ``signature`` are just standard reserved properties.

        :return:
        """
        return self.___mark___

    @mark.setter
    def mark(self, mark):
        self.___mark___ = mark

    @property
    def stamp(self):
        """
        ``name``, ``mark``, ``stamp`` ``signature`` are just standard reserved properties.

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
        if not hasattr(self.___obj___, '___parameters___'):
            raise ParametersError(
                    " <Parameters> : %r has no ___parameters___ property" % self.___obj___)
        else:
            return self.___obj___.___parameters___

    @property
    def statistic(self):
        """
        statistic of an object returns a dict contain some statistic info of the object.

        It can not be used to recover an object, neither does not represent an object uniquely.
        """
        if not hasattr(self.___obj___, '___statistic___'):
            raise StatisticError(
                " <Statistic> : %r has no ___statistic___ property" % self.___obj___)
        else:
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
