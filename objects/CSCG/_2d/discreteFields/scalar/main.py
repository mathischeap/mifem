# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 8/30/2022 11:49 PM
"""
import sys

if './' not in sys.path:
    sys.path.append('./')
import numpy as np

from objects.CSCG._2d.discreteFields.base.main import _2dCSCG_DiscreteField
from objects.CSCG._2d.discreteFields.scalar.visualize.main import _2cCSCG_DS_Visualize
from objects.CSCG._2d.discreteFields.scalar.portion import _2dCSCG_DF_ScalarPortion

class _2dCSCG_DF_Scalar(_2dCSCG_DiscreteField):
    """Region wise scalar data."""

    def __init__(self, mesh, coordinates, values, name, structured=False, grid=None):
        """"""
        super(_2dCSCG_DF_Scalar, self).__init__(mesh, coordinates, values, name,
                                                structured=structured, grid=grid)
        assert self.vdim == 1, f"vdim must be 1, now it is {self.vdim}."
        self._visualize_ = _2cCSCG_DS_Visualize(self)
        self._portion_ = _2dCSCG_DF_ScalarPortion(self)


if __name__ == '__main__':
    # mpiexec -n 4 python objects/CSCG/_2d/discreteFields/scalar/main.py
    from objects.CSCG._2d.master import MeshGenerator, SpaceInvoker, FormCaller, ExactSolutionSelector

    mesh = MeshGenerator('crazy', c=0.1)([5, 5])
    # mesh = MeshGenerator('chp1',)([2,2])
    space = SpaceInvoker('polynomials')([('Lobatto', 3), ('Lobatto', 3)])
    FC = FormCaller(mesh, space)
    ES = ExactSolutionSelector(mesh)('sL:sincos1')
    f0 = FC('0-f-o', hybrid=True)

    f0.CF = ES.potential
    f0.CF.current_time = 0
    f0.discretize()

    x = np.linspace(-1,1,10)

    ds = f0.reconstruct.discrete_scalar([x, x])

    ds.visualize()
