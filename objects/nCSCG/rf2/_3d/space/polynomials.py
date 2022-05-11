# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 11:22 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from objects.nCSCG.rf2._3d.space.base import _3nCSCG_SpaceBase
from objects.CSCG._3d.spaces.polynomials import _3dCSCG_PolynomialSpace



class _3nCSCG_PolynomialSpace(_3nCSCG_SpaceBase):
    """"""

    def __init__(self, inputs):
        """

        Parameters
        ----------
        inputs :
            Can be one of following format:
                - [a, b, c] : a, b, c are positive integers. It means [('Lobatto', a), ('Lobatto', b), ('Lobatto', c)]
                - [('Lobatto', a), ('Lobatto', b), ('Lobatto', c)]

        """
        super(_3nCSCG_PolynomialSpace, self).__init__()
        self._freeze_self_()
        self._base_space_ = _3dCSCG_PolynomialSpace(inputs, None) # the default spaces for all cells.




if __name__ == '__main__':
    # mpiexec -n 4 python objects/nCSCG/rf2/_3d/space/polynomials.py
    from objects.nCSCG.rf2._3d.__tests__.Random.mesh import random_mesh_of_elements_around as rm3
    mesh = rm3(100, refinement_intensity=0.5)
    space = _3nCSCG_PolynomialSpace([3, 3, 3])
    space.mesh = mesh
