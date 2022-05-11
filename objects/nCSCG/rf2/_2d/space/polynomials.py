# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 11:10 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from objects.nCSCG.rf2._2d.space.base import _2nCSCG_SpaceBase
from objects.CSCG._2d.spaces.polynomials import _2dCSCG_PolynomialSpace



class _2nCSCG_PolynomialSpace(_2nCSCG_SpaceBase):
    """"""

    def __init__(self, inputs):
        """

        Parameters
        ----------
        inputs :
            Can be one of following format:
                - [a, b] : a, b are two positive integers. It means [('Lobatto', a), ('Lobatto', b)]
                - [('Lobatto', a), ('Lobatto', b)]

        """
        super(_2nCSCG_PolynomialSpace, self).__init__()
        self._freeze_self_()
        self._base_space_ = _2dCSCG_PolynomialSpace(inputs, None) # the default spaces for all cells.













if __name__ == '__main__':
    # mpiexec -n 4 python objects/nCSCG/rf2/_2d/space/polynomials.py
    from objects.nCSCG.rf2._2d.__tests__.Random.mesh import random_mesh_of_elements_around as rm2
    mesh = rm2(100, refinement_intensity=0.5)
    space = _2nCSCG_PolynomialSpace([3, 3])
    space.mesh = mesh
