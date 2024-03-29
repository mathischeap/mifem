# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
import sys
from abc import ABC

if './' not in sys.path:
    sys.path.append('./')

from objects.CSCG._3d.forms.standard._1s.discretize.main import _3dCSCG_Discretize
from objects.CSCG._3d.forms.standard.base.main import _3dCSCG_Standard_Form
from objects.CSCG._3d.forms.standard._1s.project.main import _1Form_Projection
from objects.CSCG._3d.forms.standard._1s.special.main import _1Form_Special
from objects.CSCG._3d.forms.standard._1s.reconstruct import _3dCSCG_SF1_reconstruct
from objects.CSCG._3d.forms.standard._1s.inheriting.private import _3dCSCG_S1F_Private
from objects.CSCG._3d.forms.standard._1s.visualize.main import _3dCSCG_S1F_VISUALIZE
from objects.CSCG._3d.forms.standard._1s.boundary_integration.main import _3dCSCG_S1F_BI
from objects.CSCG._3d.forms.standard._1s.do.main import _3dCSCG_S1F_Do


class _3dCSCG_1Form(_3dCSCG_S1F_Private, _3dCSCG_Standard_Form, ABC):
    """
    Standard 1-form.

    :param mesh:
    :param space:
    :param hybrid:
    :param orientation:
    :param numbering_parameters:
    :param name:
    """
    def __init__(
        self, mesh, space, hybrid=True,
        orientation='outer', numbering_parameters='Naive', name=None
    ):
        if name is None:
            if hybrid:
                name = 'hybrid-' + orientation + '-oriented-1-form'
            else:
                name = orientation + '-oriented-1-form'
        super().__init__(mesh, space, hybrid, orientation, numbering_parameters, name)
        self._k_ = 1
        self.standard_properties.___PRIVATE_add_tag___('3dCSCG_standard_1form')
        self._special_ = _1Form_Special(self)
        self._projection_ = _1Form_Projection(self)
        self._discretize_ = None
        self._reconstruct_ = None
        self._visualize_ = None
        self.__BI__ = None
        self._DO_ = _3dCSCG_S1F_Do(self)
        self.___kwargs___ = {
            'hybrid': hybrid,
            'orientation': orientation,
            'numbering_parameters': numbering_parameters,
            'name': name,
        }
        self._freeze_self_()

    def ___Pr_check_CF___(self, func_body):
        assert func_body.mesh.domain == self.mesh.domain
        assert func_body.ndim == self.ndim == 3

        if func_body.__class__.__name__ == '_3dCSCG_VectorField':
            assert func_body.ftype in ('standard',), \
                f"3dCSCG 1form FUNC does not accept func _3dCSCG_VectorField of ftype {func_body.ftype}."
        else:
            raise Exception(f"3dCSCG 1form FUNC does not accept func {func_body.__class__}")

    def ___Pr_check_BC_CF___(self, func_body):
        assert func_body.mesh.domain == self.mesh.domain
        assert func_body.ndim == self.ndim == 3

        if func_body.__class__.__name__ == '_3dCSCG_VectorField':
            assert func_body.ftype in ('standard',), \
                f"3dCSCG 1form BC does not accept func _3dCSCG_VectorField of ftype {func_body.ftype}."
        else:
            raise Exception(f"3dCSCG 1form BC does not accept func {func_body.__class__}")

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
            self._reconstruct_ = _3dCSCG_SF1_reconstruct(self)
        return self._reconstruct_

    @property
    def visualize(self):
        if self._visualize_ is None:
            self._visualize_ = _3dCSCG_S1F_VISUALIZE(self)
        return self._visualize_

    @property
    def _BI_(self):
        """This properties should be accessed through `self.do.boundary_integrate`."""
        if self.__BI__ is None:
            self.__BI__ = _3dCSCG_S1F_BI(self)
        return self.__BI__


if __name__ == '__main__':
    # mpiexec -n 4 python objects/CSCG/_3d/forms/standard/_1s/main.py
    from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller

    mesh = MeshGenerator('cuboid_periodic', region_layout=[2, 2, 2])([2, 2, 2])
    space = SpaceInvoker('polynomials')([1, 1, 1])
    FC = FormCaller(mesh, space)

    # def u(t,x,y,z): return np.sin(np.pi*x)*np.cos(2*np.pi*y)*np.cos(np.pi*z) + t
    # def v(t,x,y,z): return np.cos(np.pi*x)*np.sin(np.pi*y)*np.cos(2*np.pi*z) + t
    # def w(t,x,y,z): return np.cos(np.pi*x)*np.cos(np.pi*y)*np.sin(2*np.pi*z) + t
    #
    # velocity = FC('vector', (u,v,w))
    # U = FC('scalar', u)
    # V = FC('scalar', v)
    # W = FC('scalar', w)

    f1 = FC('1-f', is_hybrid=False)
    f2 = FC('2-f', is_hybrid=False)

    E21 = f1.matrices.incidence
    E21.gathering_matrices = (f2, f1)

    E21.customize.identify_global_row(0)
    E21.customize.identify_global_row(1)
    E21.customize.identify_global_row(2)
    # E21a = E21.assembled
    print(E21.condition.condition_number)
