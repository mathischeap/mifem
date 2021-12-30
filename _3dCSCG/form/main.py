# -*- coding: utf-8 -*-
"""

A BASE for all forms except continuous forms.

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
import sys
if './' not in sys.path: sys.path.append('./')

import numpy as np
from BASE.CSCG.form.main_BASE import CSCG_FORM_BASE


# noinspection PyAbstractClass
class _3dCSCG_FORM_BASE(CSCG_FORM_BASE):
    """
    This a parent for all forms. It initialize some fundamental properties,
    like, ``mesh``, ``space``, ``ndim``, ``p`` (basis function degree) and
    ``defaultQuadDegree (dqp)``.
    """
    def __init_subclass__(cls, ndim=3):
        super().__init_subclass__(ndim=ndim)
        cls.___ndim___ = ndim

    def __init__(self, mesh, space):
        assert mesh.__class__.__name__ == '_3dCSCG_Mesh', "Need a 3dCSCG mesh."
        assert '3dCSCG|structured|space' in space.standard_properties.stamp, "Need a 3dCSCG space."
        assert mesh.ndim == space.ndim == 3
        super().__init__(mesh, space)
        self.standard_properties.___PRIVATE_add_tag___('3dCSCG_form')
        assert self.ndim == 3, "CHECK ndim"


    @property
    def dqp(self):
        """Return the Default Quadrature degree (P) for variant integrations."""
        if self.space.__class__.__name__ == '_3dCSCG_PolynomialSpace':
            return self.p
        else:
            raise NotImplementedError()





    def __neg__(self):
        """"""

    def __sub__(self, other):
        """"""
        assert other.__class__.__name__ == self.__class__.__name__, f"forms do not match."
        assert other.mesh == self.mesh, f"meshes do not match."
        assert other.space == self.space, f"spaces do not match."

    def __add__(self, other):
        """"""
        assert other.__class__.__name__ == self.__class__.__name__, f"forms do not match."
        assert other.mesh == self.mesh, f"meshes do not match."
        assert other.space == self.space, f"spaces do not match."




if __name__ == '__main__':
    # mpiexec -n 4 python _3dCSCG\form\main.py
    from _3dCSCG.main import MeshGenerator, SpaceInvoker, FormCaller#, ExactSolutionSelector

    mesh = MeshGenerator('crazy', c=0.0)([8,8,8])
    space = SpaceInvoker('polynomials')([('Lobatto',3), ('Lobatto',3), ('Lobatto',3)])
    FC = FormCaller(mesh, space)


    def p(t, x, y, z): return - 6 * np.pi * np.sin(2 * np.pi * x) * np.sin(2 * np.pi * y) * np.sin(2 * np.pi * z) + 0 * t
    def u(t,x,y,z): return np.sin(np.pi*x)*np.cos(2*np.pi*y)*np.cos(np.pi*z) + t
    def v(t,x,y,z): return np.cos(np.pi*x)*np.sin(np.pi*y)*np.cos(2*np.pi*z) + t
    def w(t,x,y,z): return np.cos(np.pi*x)*np.cos(np.pi*y)*np.sin(2*np.pi*z) + t
    velocity = FC('vector', (u,v,w))
    scalar = FC('scalar', p)

    f1 = FC('1-f', is_hybrid=False)
    F1 = FC('1-f', is_hybrid=False)


    ffa1 = f1 + F1
    ffs1 = f1 - F1

