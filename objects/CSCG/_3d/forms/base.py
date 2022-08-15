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
from functools import lru_cache
from objects.CSCG.base.forms.base.main import CSCG_FORM_BASE


# noinspection PyAbstractClass
class _3dCSCG_FORM_BASE(CSCG_FORM_BASE):
    """
    This a parent for all forms. It initializes some fundamental properties,
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
        raise NotImplementedError()

    def __sub__(self, other):
        """"""
        assert other.__class__.__name__ == self.__class__.__name__, f"forms do not match."
        assert other.mesh == self.mesh, f"meshes do not match."
        assert other.space == self.space, f"spaces do not match."

        assert self.cochain.local is not None,  f"a of (a-b) has no cochain.local, cannot perform add operator."
        assert other.cochain.local is not None, f"b of (a-b) has no cochain.local, cannot perform add operator."

        kwargs_A = self.___define_parameters___['kwargs']
        kwargs_B = other.___define_parameters___['kwargs']
        if 'name' in kwargs_A:
            name_A = kwargs_A['name']
        else:
            name_A = 'form_A'
        if 'name' in kwargs_B:
            name_B = kwargs_B['name']
        else:
            name_B = 'form_B'

        name = name_A + "-" + name_B
        kwargs_A['name'] = name

        # noinspection PyArgumentList
        result_form = self.__class__(self.mesh, self.space, **kwargs_A)

        # MUST do: add the cochain local ------------------
        COCHAIN_LOCAL = dict()
        for e in self.mesh.elements:
            clA = self.cochain.local[e]
            clB = other.cochain.local[e]
            COCHAIN_LOCAL[e] = clA - clB

        # optional 1: add the time_wise_function -----------------
        # noinspection PyBroadException
        try:
            FB = self.TW.func.body - other.TW.func.body
            result_form.TW.func.body = FB
            result_form.TW.current_time = self.TW.current_time
            result_form.TW.do.push_all_to_instant()
        except:
            pass
        #=========================================================
        result_form.cochain.local = COCHAIN_LOCAL
        return result_form

    def __add__(self, other):
        """"""
        assert other.__class__.__name__ == self.__class__.__name__, f"forms do not match."
        assert other.mesh == self.mesh, f"meshes do not match."
        assert other.space == self.space, f"spaces do not match."

        assert self.cochain.local is not None,  f"a of (a+b) has no cochain.local, cannot perform add operator."
        assert other.cochain.local is not None, f"b of (a+b) has no cochain.local, cannot perform add operator."

        kwargs_A = self.___define_parameters___['kwargs']
        kwargs_B = other.___define_parameters___['kwargs']
        if 'name' in kwargs_A:
            name_A = kwargs_A['name']
        else:
            name_A = 'form_A'
        if 'name' in kwargs_B:
            name_B = kwargs_B['name']
        else:
            name_B = 'form_B'

        name = name_A + "+" + name_B
        kwargs_A['name'] = name

        # noinspection PyArgumentList
        result_form = self.__class__(self.mesh, self.space, **kwargs_A)

        # MUST do: add the cochain local ------------------
        COCHAIN_LOCAL = dict()
        for e in self.mesh.elements:
            clA = self.cochain.local[e]
            clB = other.cochain.local[e]
            COCHAIN_LOCAL[e] = clA + clB

        # optional 1: add the time_wise_function -----------------
        # noinspection PyBroadException
        try:
            FB = self.TW.func.body + other.TW.func.body
            result_form.TW.func.body = FB
            result_form.TW.current_time = self.TW.current_time
            result_form.TW.do.push_all_to_instant()
        except:
            pass
        #=========================================================
        result_form.cochain.local = COCHAIN_LOCAL
        return result_form


    @lru_cache(maxsize=24)
    def ___PRIVATE_element_grid_data_generator_1___(self, i, density, zoom=1):
        """We generate the data for plotting all dofs of standard forms in a mesh element #`i`."""

        mesh = self.mesh
        space = self.space
        DATA = dict()

        assert i in mesh.elements, f"mesh element #{i} is not in this core."
        nodes = space.nodes
        XL, YL, ZL = nodes
        Bx = len(XL) - 1
        By = len(YL) - 1
        Bz = len(ZL) - 1

        D0, D1, D2 = None, None, None
        D0b, D1b, D2b = None, None, None
        for I, y in enumerate(YL):
            for J, z in enumerate(ZL):
                d0, d1, d2 = self.___PRIVATE_generate_line_data____('x', (y,z), density, zoom=zoom)
                if (I == 0 or I == By) and (J == 0 or J == Bz):
                    if D0b is None:
                        D0b = d0
                        D1b = d1
                        D2b = d2
                    else:
                        D0b = np.vstack((D0b, d0))
                        D1b = np.vstack((D1b, d1))
                        D2b = np.vstack((D2b, d2))
                else:
                    if D0 is None:
                        D0 = d0
                        D1 = d1
                        D2 = d2
                    else:
                        D0 = np.vstack((D0, d0))
                        D1 = np.vstack((D1, d1))
                        D2 = np.vstack((D2, d2))

        X, Y, Z = mesh.elements[i].coordinate_transformation.mapping(D0b, D1b, D2b)
        DATA['xLines_x_B'] = X # x coordinate of x-grid lines on element boundary
        DATA['xLines_y_B'] = Y # y coordinate of x-grid lines on element boundary
        DATA['xLines_z_B'] = Z # z coordinate of x-grid lines on element boundary

        if D0 is not None:
            X, Y, Z = mesh.elements[i].coordinate_transformation.mapping(D0, D1, D2)
            DATA['xLines_x'] = X # x coordinate of internal x-grid lines
            DATA['xLines_y'] = Y # y coordinate of internal x-grid lines
            DATA['xLines_z'] = Z # z coordinate of internal x-grid lines

        D0, D1, D2 = None, None, None
        D0b, D1b, D2b = None, None, None
        for I, z in enumerate(ZL):
            for J, x in enumerate(XL):
                d0, d1, d2 = self.___PRIVATE_generate_line_data____('y', (z,x), density, zoom=zoom)
                if (I == 0 or I == Bz) and (J == 0 or J == Bx):
                    if D0b is None:
                        D0b = d0
                        D1b = d1
                        D2b = d2
                    else:
                        D0b = np.vstack((D0b, d0))
                        D1b = np.vstack((D1b, d1))
                        D2b = np.vstack((D2b, d2))
                else:
                    if D0 is None:
                        D0 = d0
                        D1 = d1
                        D2 = d2
                    else:
                        D0 = np.vstack((D0, d0))
                        D1 = np.vstack((D1, d1))
                        D2 = np.vstack((D2, d2))

        X, Y, Z = mesh.elements[i].coordinate_transformation.mapping(D0b, D1b, D2b)
        DATA['yLines_x_B'] = X
        DATA['yLines_y_B'] = Y
        DATA['yLines_z_B'] = Z
        if D0 is not None:
            X, Y, Z = mesh.elements[i].coordinate_transformation.mapping(D0, D1, D2)
            DATA['yLines_x'] = X
            DATA['yLines_y'] = Y
            DATA['yLines_z'] = Z

        D0, D1, D2 = None, None, None
        D0b, D1b, D2b = None, None, None
        for I, x in enumerate(XL):
            for J, y in enumerate(YL):
                d0, d1, d2 = self.___PRIVATE_generate_line_data____('z', (x,y), density, zoom=zoom)
                if (I == 0 or I == Bx) and (J == 0 or J == By):
                    if D0b is None:
                        D0b = d0
                        D1b = d1
                        D2b = d2
                    else:
                        D0b = np.vstack((D0b, d0))
                        D1b = np.vstack((D1b, d1))
                        D2b = np.vstack((D2b, d2))
                else:
                    if D0 is None:
                        D0 = d0
                        D1 = d1
                        D2 = d2
                    else:
                        D0 = np.vstack((D0, d0))
                        D1 = np.vstack((D1, d1))
                        D2 = np.vstack((D2, d2))

        X, Y, Z = mesh.elements[i].coordinate_transformation.mapping(D0b, D1b, D2b)
        DATA['zLines_x_B'] = X
        DATA['zLines_y_B'] = Y
        DATA['zLines_z_B'] = Z

        if D0 is not None:
            X, Y, Z = mesh.elements[i].coordinate_transformation.mapping(D0, D1, D2)
            DATA['zLines_x'] = X
            DATA['zLines_y'] = Y
            DATA['zLines_z'] = Z

        #-------- below, we find the center -------------

        xyz = np.array([0]), np.array([0]), np.array([0])
        xyz = mesh.elements[i].coordinate_transformation.mapping(*xyz)
        DATA['center'] = {'coordinate': xyz, 'number': 'ME-'+str(i)}

        return DATA

    @staticmethod
    def ___PRIVATE_generate_line_data____(direction, coordinates, density, zoom=1):
        """

        :param direction:
        :param density:
        :return:
        """
        ___ = np.linspace(-1,1,density) * zoom
        if direction == 'x':
            y, z = coordinates
            X = ___
            Y = y * np.ones(density) * zoom
            Z = z * np.ones(density) * zoom
        elif direction == 'y':
            z, x = coordinates
            X = x * np.ones(density) * zoom
            Y = ___
            Z = z * np.ones(density) * zoom
        elif direction == 'z':
            x, y = coordinates
            X = x * np.ones(density) * zoom
            Y = y * np.ones(density) * zoom
            Z = ___
        else:
            raise Exception()

        return X, Y, Z








if __name__ == '__main__':
    # mpiexec -n 4 python _3dCSCG\forms\base.py
    from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller#, ExactSolutionSelector

    mesh = MeshGenerator('crazy', c=0.0)([5,5,5])
    space = SpaceInvoker('polynomials')([('Lobatto',3), ('Lobatto',3), ('Lobatto',3)])
    FC = FormCaller(mesh, space)

    def p(t, x, y, z): return - 6 * np.pi * np.sin(2 * np.pi * x) * np.sin(2 * np.pi * y) * np.sin(2 * np.pi * z) + 0 * t
    def u(t,x,y,z): return np.sin(np.pi*x)*np.cos(2*np.pi*y)*np.cos(np.pi*z) + t
    def v(t,x,y,z): return np.cos(np.pi*x)*np.sin(np.pi*y)*np.cos(2*np.pi*z) + t
    def w(t,x,y,z): return np.cos(np.pi*x)*np.cos(np.pi*y)*np.sin(2*np.pi*z) + t
    velocity = FC('vector', (u,v,w))
    scalar = FC('scalar', p)

    f1 = FC('1-f', is_hybrid=False, name='A')
    F1 = FC('1-f', is_hybrid=False)

    f1.TW.func.do.set_func_body_as(velocity)
    f1.TW.current_time = 0
    f1.TW.___DO_push_all_to_instant___()
    f1.discretize()
    F1.TW.func.do.set_func_body_as(velocity)
    F1.TW.current_time = 0
    F1.TW.___DO_push_all_to_instant___()
    F1.discretize()

    ffa1 = f1 + F1
    ffs1 = f1 - F1

    print(ffa1.error.L())
    print(ffs1.error.L())