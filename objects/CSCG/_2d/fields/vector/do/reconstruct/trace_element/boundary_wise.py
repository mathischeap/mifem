# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/08/29 7:49 PM
"""
import sys
if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
import numpy as np


class OnTraceElement_BoundaryWise(FrozenOnly):
    def __init__(self, vf):
        self._vf_ = vf
        self._freeze_self_()

    def __call__(self, xi, eta, ravel, i):
        """When on trace-element reconstruction, we only accept 1d xi, eta.

        :param xi: 1d array
        :param eta: 1d array
        :param ravel:
        :param i:
            1) self.ftype == 'standard':
                Do the reconstruction in mesh element #i. When it is None, it means all local mesh
                elements.
        :return:
        """

        SELF = self._vf_

        xyz = dict()
        value = dict()

        func = SELF.___DO_evaluate_func_at_time___()

        if i in (None, 'on_mesh_boundaries'):
            INDICES = list()
            RTE = SELF.mesh.boundaries.range_of_trace_elements
            for bn in SELF.func:  # this may not contain all mesh boundaries, only valid ones.
                INDICES.extend(RTE[bn])

        else:
            raise NotImplementedError(f"_3dCSCG_VectorField of 'boundary-wise' ftype"
                                      f" trace-element-reconstruction currently don't accept i={i}.")

        for I in INDICES:

            te = SELF.mesh.trace.elements[I]
            assert te.IS.on_mesh_boundary, f"must be the case because ftype == 'boundary-wise!"


            if te.normal_direction == 'UD':
                # only accept 1d arrays
                xyz_i = te.coordinate_transformation.mapping(xi)
            else:
                xyz_i = te.coordinate_transformation.mapping(eta)

            bn = te.on_mesh_boundary
            assert bn in func, f"trace element #{I} is on <{bn}> which is not covered by boundary-wise func."
            func_i = func[bn]

            vx_i = func_i[0](*xyz_i)
            vy_i = func_i[1](*xyz_i)

            if ravel:
                xyz[I] = [_.ravel('F') for _ in xyz_i]
                value[I] = [vx_i.ravel('F'), vy_i.ravel('F')]
            else:
                xyz[I] = xyz_i
                value[I] = [vx_i, vy_i, ]

        return xyz, value

if __name__ == "__main__":
    # mpiexec -n 4 python objects/CSCG/_2d/fields/vector/do/reconstruct/trace_element/boundary_wise.py
    from objects.CSCG._2d.master import MeshGenerator, SpaceInvoker, FormCaller#, ExactSolutionSelector

    mesh = MeshGenerator('crazy', c=0.)([3, 3])
    space = SpaceInvoker('polynomials')([1,1])
    FC = FormCaller(mesh, space)

    BV = {'Upper': [0, 0],
          'Down': [0, 0],
          'Left': [0, 0],
          'Right': [1, 0],}

    V = FC('vector', BV, name='boundary-vector')
    V.current_time = 0

    xi = np.linspace(-1,1,5)
    et = np.linspace(-1,1,4)

    xy, R = V.do.reconstruct(xi, et)

    for i in xy:
        print(R[i])