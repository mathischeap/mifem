# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
import sys
if './' not in sys.path: sys.path.append('./')
from abc import ABC

import numpy as np

from objects.CSCG._3d.forms.edge.base.main import _3dCSCG_Edge
from objects.CSCG._3d.forms.edge._0eg.discretize.main import _3dCSCG_Edge0Form_Discretize


class _3dCSCG_0Edge(_3dCSCG_Edge, ABC):
    """
    Edge 0-form.

    :param mesh:
    :param space:
    :param orientation:
    :param numbering_parameters:
    :param name:
    """
    def __init__(self, mesh, space, orientation='outer',
        numbering_parameters='Naive', name='outer-oriented-0-edge-form'):
        super().__init__(mesh, space, orientation, numbering_parameters, name)
        self._k_ = 0
        self.standard_properties.___PRIVATE_add_tag___('3dCSCG_edge_0form')
        self.RESET_cache()
        self._discretize_ = _3dCSCG_Edge0Form_Discretize(self)
        self._freeze_self_()


    def RESET_cache(self):
        super().RESET_cache()


    def ___Pr_check_CF___(self, func_body):
        assert func_body.mesh.domain == self.mesh.domain
        assert func_body.ndim == self.ndim == 3

        if func_body.__class__.__name__ == '_3dCSCG_ScalarField':
            assert func_body.ftype in ('standard',), \
                f"3dCSCG 0edge FUNC do not accept func _3dCSCG_ScalarField of ftype {func_body.ftype}."
        else:
            raise Exception(f"3dCSCG 0form FUNC do not accept func {func_body.__class__}")

    def ___Pr_check_BC_CF___(self, func_body):
        assert func_body.mesh.domain == self.mesh.domain
        assert func_body.ndim == self.ndim == 3

        if func_body.__class__.__name__ == '_3dCSCG_ScalarField':
            assert func_body.ftype in ('standard','boundary-wise'), \
                f"3dCSCG 0edge BC do not accept func _3dCSCG_ScalarField of ftype {func_body.ftype}."
        else:
            raise Exception(f"3dCSCG 1edge BC do not accept func {func_body.__class__}")


    @property
    def discretize(self):
        return self._discretize_


    def reconstruct(self, xi, eta, sigma, i=None):
        """

        :param xi: 1d array in [-1, 1]
        :param eta:  1d array in [-1, 1]
        :param sigma:  1d array in [-1, 1]
        :param i: Reconstruct in edge-element #i; when it is None, reconstruct in all edge elements.
        :return: two dictionaries.
        """

        if i is None:
            indices = self.mesh.edge.elements._locations_.keys()
        else:
            if not isinstance(i, (list, tuple)):
                indices = [i, ]
            else:
                indices = i

        basis = self.do.evaluate_basis_at_meshgrid(xi, eta, sigma)
        xyz = dict()
        v = dict()

        for i in indices:
            if i in self.mesh.edge.elements:
                ee = self.mesh.edge.elements[i]
                mesh_element = ee.CHARACTERISTIC_element
                corner_edge = ee.CHARACTERISTIC_corner_edge
                xyz_i = ee.coordinate_transformation.mapping(xi, eta, sigma,
                                                             from_element=mesh_element,
                                                             corner_edge=corner_edge)
                prime_cochain = self.cochain.local_EEW[i]
                vi = np.einsum('i, ij -> j', prime_cochain,
                               basis[corner_edge][0], optimize='greedy')
                xyz[i] = xyz_i
                v[i] = [vi, ]

        return xyz, v







if __name__ == '__main__':
    # mpiexec -n 6 python _3dCSCG\forms\edge\_0_edge.py
    from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller#, ExactSolutionSelector

    mesh = MeshGenerator('crazy', c=0.25)([5,6,7])
    space = SpaceInvoker('polynomials')([('Lobatto',10), ('Lobatto',10), ('Lobatto',10)])
    FC = FormCaller(mesh, space)

    e0 = FC('0-e')

    def p(t, x, y, z): return - 6 * np.pi * np.sin(2*np.pi*x) * np.sin(2*np.pi*y) * np.sin(2*np.pi*z) + 0 * t
    scalar = FC('scalar', p)

    e0.TW.func.do.set_func_body_as(scalar)
    e0.TW.current_time = 0
    e0.TW.do.push_all_to_instant()

    e0.discretize()
    print(e0.error.L())