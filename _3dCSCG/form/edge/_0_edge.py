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

from _3dCSCG.form.edge.base.main import _3dCSCG_Edge



class _0Edge(_3dCSCG_Edge, ABC):
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
        self._freeze_self_()


    def RESET_cache(self):
        super().RESET_cache()


    def ___TW_FUNC_body_checker___(self, func_body):
        assert func_body.mesh.domain == self.mesh.domain
        assert func_body.ndim == self.ndim == 3

        if func_body.__class__.__name__ == '_3dCSCG_ScalarField':
            assert func_body.ftype in ('standard',), \
                f"3dCSCG 0edge FUNC do not accept func _3dCSCG_ScalarField of ftype {func_body.ftype}."
        else:
            raise Exception(f"3dCSCG 0form FUNC do not accept func {func_body.__class__}")

    def ___TW_BC_body_checker___(self, func_body):
        assert func_body.mesh.domain == self.mesh.domain
        assert func_body.ndim == self.ndim == 3

        if func_body.__class__.__name__ == '_3dCSCG_ScalarField':
            assert func_body.ftype in ('standard','boundary-wise'), \
                f"3dCSCG 0edge BC do not accept func _3dCSCG_ScalarField of ftype {func_body.ftype}."
        else:
            raise Exception(f"3dCSCG 1edge BC do not accept func {func_body.__class__}")




    def discretize(self, update_cochain=True, target='func'):
        """
        Discretize the current function (a scalar field) to cochain.

        It is actually a wrapper of multiple methods that discretize functions of different types (a scalar
        field can be defined and represented in different ways in `python`, right?).

        :param bool update_cochain: (`default`: ``True``) If we update cochain with the output? Sometimes we
            may do not want to do so since we just want to use this method do some external jobs.
        :param target:
        :return: The cochain.
        :rtype: Its type can be different according to the particular discretize method.
        """
        if target == 'func':
            if self.TW.func.body.__class__.__name__ == '_3dCSCG_ScalarField':
                if self.func.ftype == 'standard':
                    return self.___PRIVATE_discretize_standard_ftype___(update_cochain=update_cochain)
                else:
                    raise NotImplementedError(f"3dCSCG 0-edge cannot (target func) discretize _3dCSCG_ScalarField of ftype={self.func.ftype}")

            else:
                raise NotImplementedError(f'3dCSCG 0-edge can not (target func) discretize {self.TW.func.body.__class__}.')

        else:
            raise NotImplementedError(f"3dCSCG 0-edge cannot discretize while targeting at {target}.")



    def ___PRIVATE_discretize_standard_ftype___(self, update_cochain=True, target='func'):
        """Discretize the standard _3dCSCG_ScalarField to a 0-edge-form

        'locally full local EEW cochain' means the cochain is a dict whose keys are edge-element
        numbers and values are edge-element-wise local cochains.

        """
        nodes = self.space.nodes

        local_EEW = dict()
        if target == 'func':
            FUNC = self.func.body[0]
        elif target == 'BC':
            FUNC = self.BC.body[0]
            assert update_cochain is False, \
                f"When target is {target}, cannot update cochain!"
        else:
            raise NotImplementedError(
                f"_0Form.___PRIVATE_discretize_standard_ftype___ "
                f"does not work for target={target}.")
        for i in self.mesh.edge.elements: # go through all local edge-elements
            element = self.mesh.edge.elements[i] # edge-element
            mesh_element = element.CHARACTERISTIC_element
            corner_edge = element.CHARACTERISTIC_corner_edge
            xyz = element.coordinate_transformation.mapping(*nodes,
                                                            from_element=mesh_element,
                                                            corner_edge=corner_edge)

            local_EEW[i] = FUNC(*xyz)
        # isKronecker? ...
        if not self.space.IS_Kronecker: raise NotImplementedError()
        # pass to cochain.local ...
        if update_cochain: self.cochain.local_EEW = local_EEW
        # ...
        return 'locally full local EEW cochain', local_EEW



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

        basis = self.DO.evaluate_basis_at_meshgrid(xi, eta, sigma)
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
    # mpiexec -n 6 python _3dCSCG\form\edge\_0_edge.py
    from _3dCSCG.main import MeshGenerator, SpaceInvoker, FormCaller#, ExactSolutionSelector

    mesh = MeshGenerator('crazy', c=0.25)([5,6,7])
    space = SpaceInvoker('polynomials')([('Lobatto',10), ('Lobatto',10), ('Lobatto',10)])
    FC = FormCaller(mesh, space)

    e0 = FC('0-e')

    def p(t, x, y, z): return - 6 * np.pi * np.sin(2*np.pi*x) * np.sin(2*np.pi*y) * np.sin(2*np.pi*z) + 0 * t
    scalar = FC('scalar', p)

    e0.TW.func.DO.set_func_body_as(scalar)
    e0.TW.current_time = 0
    e0.TW.DO.push_all_to_instant()

    e0.discretize()
    print(e0.error.L())