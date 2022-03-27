# -*- coding: utf-8 -*-
import sys
if './' not in sys.path: sys.path.append('../')
import pickle
from root.config.main import rAnk, mAster_rank
from screws.freeze.base import FrozenOnly
from screws.freeze.standard_properties import StandardProperties



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
        """Better be called from `mifem.save` when save a object."""
        _2bs_ = dict()
        # _2bs_['obj'] = str(self).split()[0][1:]
        _2bs_['obj'] = self.__class__.__name__
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
        """For a particular class, we can define its saving info by overriding this method.

        When we save multiple instances in a list or tuple at once, we will attach their
        `___PRIVATE_saving_info___` at the end of the list or tuple. Currently, this information has
        no impact. This function is just reserved, for future extensions.
        """
        return None