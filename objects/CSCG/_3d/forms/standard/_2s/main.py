# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
import sys
if './' not in sys.path: sys.path.append('./')

from root.config.main import *
from objects.CSCG._3d.forms.standard._2s.discretize.main import _3dCSCG_Discretize
from objects.CSCG._3d.forms.standard.base.main import _3dCSCG_Standard_Form
from objects.CSCG._3d.forms.standard._2s.special.main import _2Form_Special
from objects.CSCG._3d.forms.standard._2s.project.main import _2Form_Projection
from objects.CSCG._3d.forms.standard._2s.reconstruct import _3dCSCG_SF2_reconstruct
from objects.CSCG._3d.forms.standard._2s.inheriting.private import _3dCSCG_S2F_Private
from objects.CSCG._3d.forms.standard._2s.visualize.main import _3dCSCG_S2F_VISUALIZE



class _3dCSCG_2Form(_3dCSCG_S2F_Private, _3dCSCG_Standard_Form):
    """
    Standard 2-form.

    :param mesh:
    :param space:
    :param is_hybrid:
    :param orientation:
    :param numbering_parameters:
    :param name:
    """
    def __init__(self, mesh, space, is_hybrid=True,
        orientation='outer', numbering_parameters='Naive',  name=None):
        if name is None:
            if is_hybrid:
                name = 'hybrid-' + orientation + '-oriented-2-form'
            else:
                name = orientation + '-oriented-2-form'
        super().__init__(mesh, space, is_hybrid, orientation, numbering_parameters, name)
        self._k_ = 2
        self.standard_properties.___PRIVATE_add_tag___('3dCSCG_standard_2form')
        self._special_ = _2Form_Special(self)
        self._projection_ = _2Form_Projection(self)
        self.___PRIVATE_reset_cache___()
        self._discretize_ = None
        self._reconstruct_ = None
        self._visualize_ = None
        self._freeze_self_()

    def ___PRIVATE_reset_cache___(self):
        super().___PRIVATE_reset_cache___()

    def ___PRIVATE_TW_FUNC_body_checker___(self, func_body):
        assert func_body.mesh.domain == self.mesh.domain
        assert func_body.ndim == self.ndim == 3

        if func_body.__class__.__name__ == '_3dCSCG_VectorField':
            assert func_body.ftype in ('standard',), \
                f"3dCSCG 2form FUNC do not accept func _3dCSCG_VectorField of ftype {func_body.ftype}."
        else:
            raise Exception(f"3dCSCG 2form FUNC do not accept func {func_body.__class__}")

    def ___PRIVATE_TW_BC_body_checker___(self, func_body):
        assert func_body.mesh.domain == self.mesh.domain
        assert func_body.ndim == self.ndim == 3

        if func_body.__class__.__name__ == '_3dCSCG_VectorField':
            assert func_body.ftype in ('standard','boundary-wise'), \
                f"3dCSCG 2form BC do not accept func _3dCSCG_VectorField of ftype {func_body.ftype}."
        else:
            raise Exception(f"3dCSCG 2form BC do not accept func {func_body.__class__}")


    @property
    def special(self):
        return self._special_

    @property
    def projection(self):
        """A wrapper of all projection methods."""
        return self._projection_

    @property
    def discretize(self):
        if self._discretize_ is None:
            self._discretize_ = _3dCSCG_Discretize(self)
        return self._discretize_




    @property
    def reconstruct(self):
        if self._reconstruct_ is None:
            self._reconstruct_ = _3dCSCG_SF2_reconstruct(self)
        return self._reconstruct_

    @property
    def visualize(self):
        if self._visualize_ is None:
            self._visualize_ = _3dCSCG_S2F_VISUALIZE(self)
        return self._visualize_





if __name__ == '__main__':
    # mpiexec -n 5 python objects/CSCG/_3d/forms/standard/_2s/main.py
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
    f2.visualize(x=0.25)