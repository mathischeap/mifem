# -*- coding: utf-8 -*-

import sys
if './' not in sys.path: sys.path.append('./')

from objects.CSCG._3d.forms.standard.base.visualize.main import _3dCSCG_FormVisualize
from objects.CSCG._3d.forms.standard._2s.visualize.matplot import _3dCSCG_S2F_VISUALIZE_Matplot

import numpy as np

class _3dCSCG_S2F_VISUALIZE(_3dCSCG_FormVisualize):
    """"""
    def __init__(self, sf):
        """"""
        super(_3dCSCG_S2F_VISUALIZE, self).__init__(sf)
        self._sf_ = sf
        self._matplot_ = None
        self._paraview_ = None
        self._freeze_self_()

    def __call__(self, *args, **kwargs):
        return self.matplot(*args, **kwargs)

    @property
    def matplot(self):
        if self._matplot_ is None:
            self._matplot_ = _3dCSCG_S2F_VISUALIZE_Matplot(self._sf_)
        return self._matplot_


if __name__ == '__main__':
    # mpiexec -n 4 python objects/CSCG/_3d/forms/standard/_2s/visualize/main.py
    from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller

    mesh = MeshGenerator('crazy', c=0.0)([3,3,3])
    space = SpaceInvoker('polynomials')([('Lobatto',3), ('Lobatto',3), ('Lobatto',3)])
    FC = FormCaller(mesh, space)

    def u(t,x,y,z): return np.sin(np.pi*x)*np.cos(2*np.pi*y)*np.cos(np.pi*z) + t
    def v(t,x,y,z): return np.cos(np.pi*x)*np.sin(np.pi*y)*np.cos(2*np.pi*z) + t
    def w(t,x,y,z): return np.cos(np.pi*x)*np.cos(np.pi*y)*np.sin(2*np.pi*z) + t

    velocity = FC('vector', (u,v,w))
    U = FC('scalar', u)
    V = FC('scalar', v)
    W = FC('scalar', w)

    f2 = FC('2-f', is_hybrid=False)
    f2.TW.func.do.set_func_body_as(velocity)
    f2.TW.current_time = 0
    f2.TW.___DO_push_all_to_instant___()
    f2.discretize()