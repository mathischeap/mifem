# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/13 2:55 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.main import FrozenClass


class nCSCG_FieldBase(FrozenClass):
    """"""

    def __init__(self, mesh, ftype='standard', valid_time=None, name='scalar-field'):
        """"""
        self._func_ = None
        self._ftype_ = None
        self._mesh_ = mesh
        assert ftype in ("standard",), f"ftype={ftype} wrong." # allowed ftype.
        self.___PRIVATE_set_vt_to___(valid_time)
        self.standard_properties.___PRIVATE_add_tag___('nCSCG_RF2_Field')
        self.standard_properties.name = name
        self._current_time_ = None

    @property
    def func(self):
        """The function body: same in all cores!"""
        return self._func_

    @property
    def ftype(self):
        """The function type: same in all cores!"""
        return self._ftype_

    @property
    def mesh(self):
        return self._mesh_

    @property
    def ndim(self):
        return self.mesh.ndim


    #----------- VALID TIME -----------------------------------------------------------------
    @property
    def valid_time(self):
        """``current time`` must be with in ``valid_time``.

        Some continuous forms are only valid for some certain times.

        `valid_time` cannot be changed once we have defined the instance.
        """
        return self._valid_time_

    def ___PRIVATE_set_vt_to___(self, valid_time):
        """
        :param valid_time:
            - None                             : It can be everything and be changed whenever you want.
            - 'valid_only_at_its_first_instant': as it says...
            - int or float                     : Can only be this particular time instance.
        :return:
        """
        if valid_time is None:
            # can be any time.
            pass
        elif valid_time == 'valid_only_at_its_first_instant':
            # as the string says, we can only set its ``current_time`` once.
            pass
        elif isinstance(valid_time, (int, float)):
            # current_time can only be this single time instance: `current_time` = `valid_time`.
            pass
        else:
            raise Exception(f'valid_time = {valid_time} format wrong.')

        self._valid_time_ = valid_time

    def ___PRIVATE_check_ct_in_vt___(self, ct):
        if self.valid_time is None:
            # valid at any given time instant (current time).
            return
        elif self.valid_time == 'valid_only_at_its_first_instant':
            # current time is only valid at its first set instant.
            if self._current_time_ is None:
                pass
            elif ct == self._current_time_:
                pass
            else:
                raise Exception(f"time is `valid_only_at_its_first_instant`. "
                                f"current_time = {self._current_time_}, "
                                f"can not change it to {ct}.")
        elif isinstance(self.valid_time, (int, float)):
            if ct != self.valid_time:
                raise Exception(f"{self} can only be valid at time {self.valid_time}."
                                f" Now set to {ct}.")
        else:
            raise Exception(f'current_time = {ct} ({ct.__class__.__name__}) '
                            f'(valid_time={self.valid_time}) is illegal.')

    # ...........CURRENT TIME which is related to valid time .....................................
    @property
    def current_time(self):
        """The current time. If push to instant, we use this time."""
        assert self._current_time_ is not None, "current_time is None, set it firstly."
        return self._current_time_

    @current_time.setter
    def current_time(self, current_time):
        assert isinstance(current_time, (float,int))
        self.___PRIVATE_check_ct_in_vt___(current_time)
        self._current_time_ = current_time

    #====================== BELOW: children must have these methods =======================
    def ___Pr_set_func___(self, *args, **kwargs):
        raise NotImplementedError()

    def ___Pr_evaluate_func___(self, *args, **kwargs):
        raise NotImplementedError()

    @property
    def shape(self):
        raise NotImplementedError()





if __name__ == "__main__":
    # mpiexec -n 4 python objects/nCSCG/rf2/base/fields/base.py
    pass
