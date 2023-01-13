# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/10/08 7:05 PM
"""

from components.freeze.base import FrozenOnly
import random
import numpy as np


class miUsTriangle_ExactSolutionBase(FrozenOnly):
    """"""

    def __init__(self, mesh):
        """"""
        self._mesh_ = mesh

    @property
    def mesh(self):
        return self._mesh_

    # ... valid time: override it to make the exact solution only valid at certain time. ...
    @property
    def valid_time(self):
        return None

    def ___PRIVATE_generate_random_valid_time_instances___(self, amount=None):
        """We will generate some random valid time instances and put them in a 1d array. They can
        be used for purposes like self-checking.

        For valid time:
            - None                             : It can be everything and be changed whenever you want.
            - 'valid_only_at_its_first_instant': as it says...
            - int or float                     : Can only be this particular time instance.

        :param amount: {None, positive int}
            How many random time instances you want?  (will be functional only when it is applicable.)
        :return:
        """
        vt = self.valid_time
        if vt is None:
            if amount is None:
                amount = random.randint(2, 5)
            rTIs = np.random.rand(amount) * 10
        elif vt == 'valid_only_at_its_first_instant':
            rTIs = np.random.rand(1) * 10
        elif isinstance(vt, (int, float)):
            rTIs = np.array([vt, ])
        else:
            raise NotImplementedError(f"valid_time = {vt} is not understandable!")

        return rTIs

    @property
    def current_time(self):
        """Return a list of current_time of all valid properties.."""
        ct = list()
        for attr_name in self.__dict__:
            attr = getattr(self, attr_name)
            if hasattr(attr, '_current_time_'):
                if attr._current_time_ is not None:
                    ct.append(attr.current_time)
                else:
                    pass
            else:
                pass
        return ct

    @current_time.setter
    def current_time(self, ct):
        """Set current time for all valid attributes."""
        attr_names = dir(self)
        for attr_name in attr_names:
            if attr_name != 'current_time':
                attr = getattr(self, attr_name)

                if hasattr(attr, '_current_time_'):

                    attr.current_time = ct

                else:
                    pass
