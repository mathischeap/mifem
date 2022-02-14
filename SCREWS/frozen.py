# -*- coding: utf-8 -*-



import sys
if './' not in sys.path: sys.path.append('./')
import pickle
import sentry_sdk
from root.config import *
if seNtry_on: sentry_sdk.init("https://79fb951c3ea7457c8c3f47c3dfb026ce@sentry.io/1458280")

from SCREWS.exceptions import StatisticError
from SCREWS.exceptions import ParametersError
from SCREWS.exceptions import FrozenError

from abc import ABC




class FrozenOnly(ABC):
    @property
    def ___isfrozen___(self):
        try:
            return self.___ISFROZEN___
        except AttributeError:
            object.__setattr__(self, '___ISFROZEN___', False)
            return self.___ISFROZEN___

    def __setattr__(self, key, value):
        if self.___isfrozen___ and key not in dir(self):
            raise FrozenError(" <FrozenClass> : %r is frozen. CANNOT define new attributes." % self)
        object.__setattr__(self, key, value)

    def _freeze_self_(self):
        """Freeze self, can add no more new attributes. """
        try:
            # we will run method ___PreFrozenChecker___() if it exists.
            getattr(self, '___PreFrozenChecker___')()
        except AttributeError:
            pass
        self.___ISFROZEN___ = True

    def _melt_self_(self):
        """Melt self, so  we can define new attributes."""
        self.___ISFROZEN___ = False

    @property
    def _frozen_(self):
        """Return the status of the form, frozen (True) or melt (False)?"""
        return self.___isfrozen___



class FrozenClass(FrozenOnly):

    @property
    def standard_properties(self):
        """
        A wrapper of all standard properties.

        It has, for example, following sub-properties:

        - **name** -- (str) The name of the object.
        - **mark** -- (str) The mark of the object.
        - **stamp** -- (str) The stamp of the object.
        - **parameters** -- (dict) The parameters of the object. To use this, we need first define a private method
            ``___parameters___`` in the object class to collect the parameters.
        - **statistic** -- (dict) The statistic of the object. To use this, we need first define a private method
            ``___statistic___`` in the object class to collect the statistic.
        """
        if not hasattr(self, '___sp___'):
            if self._frozen_:
                self._melt_self_()
                self.___sp___ = StandardProperties(self)
                self._freeze_self_()
            else:
                self.___sp___ = StandardProperties(self)
        return self.___sp___


    def ___PRIVATE_save___(self, filename, do_save=False):
        """Better be called from mifem.save when save a object."""
        _2bs_ = dict()
        _2bs_['obj'] = str(self).split()[0][1:]
        _2bs_['parameters'] = self.standard_properties.parameters
        if rAnk == mAster_rank:
            if do_save:
                if filename[-3:] != '.mi': filename += '.mi'
                assert filename.count('.') == 1, f"filename={filename} wrong."
                with open(filename, 'wb') as output:
                    pickle.dump(_2bs_, output, pickle.HIGHEST_PROTOCOL)
                output.close()
        return _2bs_

    # noinspection PyMethodMayBeStatic
    def ___PRIVATE_saving_info___(self):
        """For a particular class, we can define its saving info by overriding this method. By doing this, we can
        use this info for special reading routine when reading multiple objects at once. For example, if we read four
        forms who have the same mesh and the same space, then we do not need to rebuild this mesh (or space) for four
        times. So we have to put information into the save file when saving these four forms."""
        return None



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