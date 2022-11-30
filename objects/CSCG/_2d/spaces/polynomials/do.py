# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 11/29/2022 11:50 AM
"""
import sys

if './' not in sys.path: sys.path.append('./')
from objects.CSCG._2d.spaces.base.do import _2dCSCG_space_do


class _2dCSCG_space_Polynomial_Do(_2dCSCG_space_do):
    """"""

    def __init__(self, space):
        """"""
        super(_2dCSCG_space_Polynomial_Do, self).__init__(space)
        self._freeze_self_()


    def refine(self, p=(1, 1)):
        """Lets say self.space.p = (a, b), and input p = (c, d), this function gives a new
        space of degree (a+c, b+d)

        Parameters
        ----------
        p

        Returns
        -------

        """
        assert len(p) == 2 and all([isinstance(_, int) for _ in p]), \
            f"p={p} wrong, must be a tuple or list of 2 integers (no need to be positive)."

        basis_1d = self.space.basises
        category = [_.category for _ in basis_1d]
        px, py = self.space.p
        px += p[0]
        py += p[0]
        new_input = (
            [category[0], px],
            [category[1], py]
        )

        return self.space.__class__(new_input, None)



if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
