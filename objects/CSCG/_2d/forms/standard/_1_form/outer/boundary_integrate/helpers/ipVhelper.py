# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/08/29 7:34 PM
"""

from components.freeze.base import FrozenOnly


class ipV_Helper(FrozenOnly):
    """"""

    def __init__(self, s1f, V, quad_degree):
        """

        Parameters
        ----------
        s1f : standard-1-form
            The self 1-form. No need to have a cochain.
        V :
            A 3dCSCG vector.
        quad_degree
        """


    def __call__(self, e):
        """"""
        raise NotImplementedError()
