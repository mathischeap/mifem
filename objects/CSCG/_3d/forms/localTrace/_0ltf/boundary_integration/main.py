# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 12/19/2022 11:36 AM
"""
import sys

if './' not in sys.path:
    sys.path.append('./')
from components.freeze.main import FrozenOnly
import numpy as np
from objects.CSCG._3d.forms.localTrace._0ltf.boundary_integration.helpers.scalar import Scalar
from tools.elementwiseCache.dataStructures.objects.columnVector.main import EWC_ColumnVector


class _3dCSCG_0LTF_BoundaryIntegration(FrozenOnly):
    """"""

    def __init__(self, ltf):
        """"""
        self._ltf_ = ltf
        self._freeze_self_()

    def __call__(self, obj, **kwargs):
        """"""

        if obj.__class__.__name__ == '_3dCSCG_ScalarField':
            return self.___Pr__3dCSCG_ScalarField___(obj, **kwargs)
        else:
            raise NotImplementedError()

    def ___Pr__3dCSCG_ScalarField___(self, S, quad_degree=None):
        """"""
        if quad_degree is None:
            quad_degree = [self._ltf_.dqp[i] + 1 for i in range(self._ltf_.ndim)]
        else:
            pass

        dg = Scalar(self._ltf_, S, quad_degree)

        return EWC_ColumnVector(self._ltf_.mesh.elements, dg, 'no_cache')


if __name__ == '__main__':
    # mpiexec -n 4 python objects/CSCG/_3d/forms/localTrace/_0ltf/boundary_integration/main.py
    from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller

    mesh = MeshGenerator('ct', c=0.)([3, 3, 3])
    space = SpaceInvoker('polynomials')([2, 2, 2])
    FC = FormCaller(mesh, space)

    ltf0 = FC('0-lt', hybrid=False)

    def p(t, x, y, z): return np.cos(2*np.pi*x) * np.cos(np.pi*y) * np.cos(3*np.pi*z) + t

    scalar = FC('scalar', p)
    scalar.current_time = 0

    # ltf0.BC.CF = scalar
    # ltf0.BC.boundaries = 'all'
    # BCM = ltf0.BC.interpret.boundary_cochain
    # M = ltf0.matrices.mass
    #
    # V = M @ BCM
    #
    # for i in V:
    #     print(i, V[i].toarray().T)

    bi = ltf0.boundary_integration(scalar)
