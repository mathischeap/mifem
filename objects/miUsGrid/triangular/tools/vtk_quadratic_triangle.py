# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 10/5/2022 10:51 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')
import numpy as np
from root.config.main import RANK, MASTER_RANK
from pyevtk.hl import unstructuredGridToVTK


def vtk_quadratic_triangle(*forms, filename='vtk_quadratic_triangle'):
    """export a series of 2d triangular forms into a single vtu file of vtk_quadratic_triangle cells."""

    mesh = forms[0].mesh
    assert mesh.__class__.__name__ == 'miUsGrid_TriangularMesh'
    for i, form in enumerate(forms):
        assert mesh.ndim == 2
        assert form.mesh == mesh, f"form[{i}]'s mesh does not match."

    fD = forms[0].export.vtk.quadratic_triangle(data_only=True)
    if RANK == MASTER_RANK:
        xyz, con, data = fD

    if len(forms) == 1:
        pass
    else:
        for f in forms[1:]:
            fD = f.export.vtk.quadratic_triangle(data_only=True)
            if RANK == MASTER_RANK:
                # noinspection PyUnboundLocalVariable
                data.update(fD[2])

    if RANK == MASTER_RANK:
        # noinspection PyUnboundLocalVariable
        num_cells = int(len(con) / 6)
        offsets = np.linspace(6, 6 * num_cells, num_cells, dtype=int)
        # noinspection PyUnboundLocalVariable
        unstructuredGridToVTK(
            filename,
            *xyz,
            connectivity = con, offsets=offsets,
            cell_types = 22 * np.ones_like(offsets),
            pointData=data
        )



if __name__ == '__main__':
    # mpiexec -n 4 python objects/miUsGrid/tools/vtk_quadratic_triangle.py
    from __init__ import miTri

    fc = miTri.call('st32', 3)


    def func(t, x, y): return np.sin(2 * np.pi * x) * np.cos(2 * np.pi * y) + t

    scalar = fc('scalar', func)
    scalar.current_time = 0

    f0 = fc('0-f-o', name='pressure')
    f0.CF = scalar
    f0.discretize()

    f2 = fc('2-f-i', name = 'potential')
    f2.CF = scalar
    f2.discretize()

    def p_func(t, x, y): return np.sin(2 * np.pi * x) * np.cos(2 * np.pi * y) + t
    def q_func(t, x, y): return np.cos(2 * np.pi * x) * np.sin(2 * np.pi * y) + t

    v = fc('vector', [p_func,q_func])
    v.current_time = 0

    f1o = fc('1-f-o', name='velocity')
    f1o.CF = v
    f1o.discretize()
    f1i = fc('1-f-i', name='vorticity')
    f1i.CF = v
    f1i.discretize()

    vtk_quadratic_triangle(f0, f1i, f1o, f2)
