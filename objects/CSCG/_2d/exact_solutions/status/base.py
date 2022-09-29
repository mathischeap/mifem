# -*- coding: utf-8 -*-
from screws.freeze.main import FrozenOnly
import random
import numpy as np


class Base(FrozenOnly):
    """A base (parent) for all exact solution classes."""
    def __init__(self, es):
        self._es_ = es
        self._mesh_ = es.mesh

    # ... valid time: override it to make the exact solution only valid at certain time. ...
    @property
    def valid_time(self):
        return None

    @property
    def mesh(self):
        """
        all exact solutions have mesh, no matter if mesh is useful (general exact solutions not use mesh).

        :return: A ``_2dCSCG_Mesh`` object.
        """
        return self._mesh_

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
            if amount is None: amount = random.randint(2, 5)
            rTIs = np.random.rand(amount) * 10
        elif vt == 'valid_only_at_its_first_instant':
            rTIs = np.random.rand(1) * 10
        elif isinstance(vt, (int, float)):
            rTIs = np.array([vt,])
        else:
            raise NotImplementedError(f"valid_time = {vt} is not understandable!")

        return rTIs