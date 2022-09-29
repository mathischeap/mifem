# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 9/2/2022 1:07 AM
"""
import sys

if './' not in sys.path: sys.path.append('./')
import numpy as np
from evtk import hl
from screws.miscellaneous.mios import mkdir, cleandir
from root.config.main import cOmm


def gridToVTK(grid, objs, path):
    """Export the objs to a vtk dataset named "filename" on a grid defined by `grid`.

    The VTK data will be region-wise: each region has its own vtk file. Loading all these
    vtk files to a visualization tool will reveal the whole domain.

    This is very good for a domain of one single region and is okay for a domain of a few regions.

    While, for a domain of many regions, it is very inconvenient. Use `unstructuredGridToVTK`
    instead.

    Parameters
    ----------
    grid : {list, dict}
        If it is a list, we will make it into a dict whose keys are the region names  of the mesh,
        and this list is the value for all keys.

        For each value, if the mesh is 3 (2) d, it is a list of three (two) 1d array.

        The mesh will be the mesh of the `objs`. So we need all objs to be  of the same mesh.
    objs
    path

    Returns
    -------
    str
        Full path to saved file.

    """
    if not isinstance(objs, (list, tuple)):
        objs = [objs,]
    else:
        pass
    for i, obj in enumerate(objs):
        assert hasattr(obj, 'mesh'),  f"{i}th obj does not have a mesh."
        if i == 0:
            mesh = obj.mesh
        else:
            assert mesh == obj.mesh, f"mesh of {i}th obj does not match that of 0th obj."

    DiscreteFields = list()
    for obj in  objs:
        df = obj.reconstruct.discrete_field(grid)
        DiscreteFields.append(df)

    df0 = DiscreteFields[0]
    regions = df0.regions # the local regions; data for these regions are stored locally in this core.

    cOmm.barrier() # make sure the folder is really for all cores.
    mkdir(path)
    cleandir(path)

    #------ variable names ------------------------------------------------------------------------1
    var_names = list()
    for i, obj in enumerate(objs):
        if hasattr(obj, 'name'):
            name = obj.name
        elif hasattr(obj, 'standard_properties'):
            name = str(obj.standard_properties.name)
        else:
            name = f'variable_{i}'

        j = 0
        old_name = name
        while name in var_names:
            j += 1
            name = old_name + '_' + str(j)

        var_names.append(name)

    cOmm.barrier() # make sure the folder is really for all cores.
    for i, rn in enumerate(regions):
        if i == 0:
            xyz = df0.coordinates[rn]

        VALs = dict()
        for j, df in enumerate(DiscreteFields):
            if df.vdim == 1:
                VALs[var_names[j]] = df.values[rn][0]
            else:
                VALs[var_names[j]] = tuple(df.values[rn])

        if len(xyz) == 2: # we are in a 2d CSCG mesh ------
            x, y = xyz
            x = x[:,:,np.newaxis]
            y = y[:,:,np.newaxis]
            z = np.zeros_like(x)

            new_VALs = dict()
            for vn in VALs:
                if VALs[vn].__class__.__name__ == 'ndarray':
                    new_VALs[vn] = VALs[vn][:,:,np.newaxis]
                else:
                    new_VALs[vn] = list()
                    for data in VALs[vn]:
                        new_VALs[vn].append(data[:,:,np.newaxis])
                    new_VALs[vn].append(np.zeros_like(new_VALs[vn][0]))
                    new_VALs[vn] = tuple(new_VALs[vn])

            hl.gridToVTK(
                f"./{path}/" + "region_" + rn[2:],
                x, y, z,
                pointData=new_VALs
            )

        elif len(xyz) == 3:

            hl.gridToVTK(
                f"./{path}/" + "region_" + rn[2:],
                *xyz,
                pointData=VALs
            )

        else:
            raise Exception(f"Trivial check!")



if __name__ == '__main__':
    # mpiexec -n 4 python objects/CSCG/tools/gridToVTK.py
    from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller

    mesh = MeshGenerator('crazy', c=0.1)([3,3,3])

    space = SpaceInvoker('polynomials')([('Lobatto',3), ('Lobatto',3), ('Lobatto',3)])
    FC = FormCaller(mesh, space)

    def u(t,x,y,z): return np.sin(np.pi*x)*np.cos(2*np.pi*y)*np.cos(np.pi*z) + t
    def v(t,x,y,z): return np.cos(np.pi*x)*np.sin(np.pi*y)*np.cos(2*np.pi*z) + t
    def w(t,x,y,z): return np.cos(np.pi*x)*np.cos(np.pi*y)*np.sin(2*np.pi*z) + t

    def p(t,x,y,z): return np.sin(2*np.pi*x)*np.sin(np.pi*y)*np.sin(2*np.pi*z) + t

    scalar = FC('scalar', p)
    velocity = FC('vector', (u,v,w))
    U = FC('scalar', u)
    V = FC('scalar', v)
    W = FC('scalar', w)

    f0 = FC('0-f', is_hybrid=False, name='pressure')
    f0.TW.func.do.set_func_body_as(scalar)
    f0.TW.current_time = 0
    f0.TW.___DO_push_all_to_instant___()
    f0.discretize()
    f1 = FC('1-f', is_hybrid=False, name='vorticity')
    f1.TW.func.do.set_func_body_as(velocity)
    f1.TW.current_time = 0
    f1.TW.___DO_push_all_to_instant___()
    f1.discretize()
    f2 = FC('2-f', is_hybrid=False, name='velocity')
    f2.TW.func.do.set_func_body_as(velocity)
    f2.TW.current_time = 0
    f2.TW.___DO_push_all_to_instant___()
    f2.discretize()
    f3 = FC('3-f', is_hybrid=False, name='total pressure')
    f3.TW.func.do.set_func_body_as(scalar)
    f3.TW.current_time = 0
    f3.TW.___DO_push_all_to_instant___()
    f3.discretize()

    grid = [np.linspace(-1,1,10), np.linspace(-1,1,10), np.linspace(-1,1,10)]

    gridToVTK(grid, [f0, f1, f2, f3], 'gridToVTK_test')


    # ----------- 2d test below ------------------------------------------------
    from objects.CSCG._2d.master import MeshGenerator, SpaceInvoker, FormCaller

    # mesh = MeshGenerator('crazy', c=0.1)([3,3])
    #
    # space = SpaceInvoker('polynomials')([3,3])
    # FC = FormCaller(mesh, space)
    #
    # def u(t,x,y): return np.sin(np.pi*x)*np.cos(2*np.pi*y)+t
    # def v(t,x,y): return np.cos(np.pi*x)*np.sin(np.pi*y)+t
    #
    # def p(t,x,y): return np.sin(2*np.pi*x)*np.sin(np.pi*y)+t
    #
    # scalar = FC('scalar', p)
    # velocity = FC('vector', (u,v))
    #
    # f0 = FC('0-f-o', is_hybrid=False, name='pressure')
    # f0.TW.func.do.set_func_body_as(scalar)
    # f0.TW.current_time = 0
    # f0.TW.___DO_push_all_to_instant___()
    # f0.discretize()
    # f1i = FC('1-f-i', is_hybrid=False, name='Velocity')
    # f1i.TW.func.do.set_func_body_as(velocity)
    # f1i.TW.current_time = 0
    # f1i.TW.___DO_push_all_to_instant___()
    # f1i.discretize()
    # f1o = FC('1-f-o', is_hybrid=False, name='vorticity')
    # f1o.TW.func.do.set_func_body_as(velocity)
    # f1o.TW.current_time = 0
    # f1o.TW.___DO_push_all_to_instant___()
    # f1o.discretize()
    # f2 = FC('2-f-o', is_hybrid=False, name='total pressure')
    # f2.TW.func.do.set_func_body_as(scalar)
    # f2.TW.current_time = 0
    # f2.TW.___DO_push_all_to_instant___()
    # f2.discretize()
    #
    # grid = [np.linspace(-1,1,50), np.linspace(-1,1,50)]
    #
    # gridToVTK(grid, [f0, f1i, f1o, f2], 'gridToVTK_2dtest')