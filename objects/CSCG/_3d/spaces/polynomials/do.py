# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 11/28/2022 10:12 PM
"""
from objects.CSCG._3d.spaces.base.do import _3dCSCG_space_do

class _3dCSCG_space_Polynomial_do(_3dCSCG_space_do):
    """"""

    def __init__(self, space):
        """"""
        super(_3dCSCG_space_Polynomial_do, self).__init__(space)
        self._freeze_self_()

    def refine(self, p =(1, 1, 1)):
        """We return a new polynomial space of higher degree.

        For example, if current space's degree is (1,2,3) and do self.do.refine(p=(2,3,4)),
        we obtain a new space of degree (1+2, 2+3, 3+4).

        Parameters
        ----------
        p

        Returns
        -------

        """
        assert len(p) == 3 and all([isinstance(_, int) for _ in p]), \
            f"p={p} wrong, must be a tuple or list of 3 integers (no need to be positive)."

        basis_1d = self._space_.basises
        category = [_.category for _ in basis_1d]
        px, py, pz = self._space_.p
        px += p[0]
        py += p[0]
        pz += p[0]
        new_input = (
            [category[0], px],
            [category[1], py],
            [category[2], pz],
        )

        return self._space_.__class__(new_input, None)