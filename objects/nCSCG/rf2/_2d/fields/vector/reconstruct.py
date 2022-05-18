# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/15 7:45 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly




class _2nCSCG_RF2_VectorReconstruct(FrozenOnly):
    """"""

    def __init__(self, cf):
        """"""
        self._cf_ = cf
        self._freeze_self_()

    def __call__(self, xi_eta, where=None, i=None, ravel=False):
        """

        Parameters
        ----------
        xi_eta
        ravel
        i
        where

        Returns
        -------

        """
        ftype = self._cf_.ftype

        # parse where when it is None --------------------------------------------
        if where is None:
            if ftype == 'standard':
                where = 'cell'
            else:
                raise NotImplementedError(f"please set default `where` for {ftype} 2nCSCG scalar field.")
        else:
            pass

        #---- parse i -----------------------------------------------------------
        if i is None:
            INDICES = self._cf_.mesh
        else:
            raise NotImplementedError(f"i={i} is not implemented.")

        if ftype == 'standard':
            if where == 'cell':
                return self.___Pr_standard_cell___(xi_eta, ravel, INDICES)
            else:
                raise NotImplementedError(f"On {where} for standard type not coded..")
        else:
            raise NotImplementedError(f"{ftype}.")

    def ___Pr_standard_cell___(self, xi_eta, ravel, INDICES):
        """

        Parameters
        ----------
        xi_eta :
        ravel :
        INDICES :
            In those cells we will reconstruct the scalar.

        Returns
        -------

        """
        func = self._cf_.___Pr_evaluate_func___()
        mesh = self._cf_.mesh

        if INDICES is mesh:
            full = True
        else:
            full = False

        assert xi_eta.___Pr_is_2nCSCG_RF2_mesh_coo___, f"I need a coo distribution object."

        xy = dict()
        value = dict()

        for i in INDICES:
            cell = mesh(i)
            xy_i = cell.coordinate_transformation.mapping(*xi_eta[i])
            v_i = [func[0](*xy_i), func[1](*xy_i)]

            if ravel:
                xy[cell.__repr__()] = [I.ravel('F') for I in xy_i]
                value[cell.__repr__()] = [I.ravel('F') for I in v_i]
            else:
                xy[cell.__repr__()] = xy_i
                value[cell.__repr__()] = v_i

        xy = self._cf_.mesh.ids('vector', xy, 2, xi_eta.distribution, full)
        value = self._cf_.mesh.ids('vector', value, 2, xi_eta.distribution, full)

        return xy, value





if __name__ == "__main__":
    # mpiexec -n 4 python objects/nCSCG/rfT2/_2d/fields/vector/reconstruct.py
    from objects.nCSCG.rf2._2d.__tests__.Random.mesh import random_mesh_of_elements_around as rm2
    from objects.nCSCG.rf2._2d.fields.vector.main import _2nCSCG_RF2_VectorField

    mesh = rm2(100)

    import numpy as np
    def p(t, x, y): return np.sin(np.pi*x) * np.sin(np.pi*y) + t
    def q(t, x, y): return np.cos(np.pi*x) * np.cos(np.pi*y) + t

    v = _2nCSCG_RF2_VectorField(mesh, (p,q))
    v.current_time = 0

    coo = mesh.coordinates.homogeneous(10, ndim=2)

    xy, v = v.reconstruct(coo)
