# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 8/28/2022 1:55 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')
import numpy as np

from objects.CSCG._2d.discreteFields.base.main import _2dCSCG_DiscreteField
from objects.CSCG._2d.discreteFields.vector.visualize.main import _2dCSCG_DV_Visualize
from objects.CSCG._2d.discreteFields.vector.portion import _2dCSCG_DF_VectorPortion

class _2dCSCG_DF_Vector(_2dCSCG_DiscreteField):
    """Region wise vector data."""

    def __init__(self, mesh, coordinates, values, name, structured=False, grid=None):
        """"""
        super(_2dCSCG_DF_Vector, self).__init__(mesh, coordinates, values, name,
                                                structured=structured, grid=grid)
        assert self.vdim == self.mesh.ndim, f"vdim must be 2, now it is {self.vdim}."
        self._visualize_ = _2dCSCG_DV_Visualize(self)
        self._portion_ = _2dCSCG_DF_VectorPortion(self)



if __name__ == '__main__':
    # mpiexec -n 4 python objects/CSCG/_2d/discrete_fields/vector/main.py
    from objects.CSCG._2d.master import MeshGenerator, SpaceInvoker, FormCaller, ExactSolutionSelector

    # mesh = MeshGenerator('crazy', c=0.3)([50,45])
    # mesh = MeshGenerator('chp1',)([2,2])
    mesh = MeshGenerator('chp2')([10,10])
    space = SpaceInvoker('polynomials')([('Lobatto',3), ('Lobatto',3)])
    FC = FormCaller(mesh, space)

    ES = ExactSolutionSelector(mesh)('sL:sincos1')

    u = FC('1-f-o', is_hybrid=True)

    u.TW.func.do.set_func_body_as(ES, 'velocity')
    u.TW.current_time = 0
    u.TW.do.push_all_to_instant()
    u.discretize()

    r = np.linspace(-1,1,20)
    s = np.linspace(-1,1,20)

    dv = u.reconstruct.discrete_vector([r, s])
    x = [-1., -0]
    y = [-1., -0]

    portion = dv.portion.xy_range(x, y)
    dv.visualize.matplot(scale=50)
    portion.visualize.matplot.quiver(scale=100)