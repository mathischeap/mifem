# -*- coding: utf-8 -*-
"""

@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/13 4:16 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly





class mpRfT2_ScalarReconstruct(FrozenOnly):
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
        where
        i
        ravel

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
        func = self._cf_.___Pr_evaluate_func___()[0]
        mesh = self._cf_.mesh

        if INDICES is mesh:
            full = True
        else:
            full = False

        assert xi_eta.___Pr_is_mpRfT2_mesh_coo_map___, f"I need a mpRfT2 coo_map object."

        xy = dict()
        value = dict()

        for i in INDICES:
            cell = mesh[i]
            xy_i = cell.coordinate_transformation.mapping(*xi_eta[i])
            v_i = func(*xy_i)

            if ravel:
                xy[cell.__repr__()] = [I.ravel('F') for I in xy_i]
                value[cell.__repr__()] = [v_i.ravel('F'),]
            else:
                xy[cell.__repr__()] = xy_i
                value[cell.__repr__()] = [v_i,]

        xy    = mesh.rcWds.vector(   xy, 2, xi_eta.distribution, full)
        value = mesh.rcWds.scalar(value, 2, xi_eta.distribution, full)

        return xy, value









if __name__ == "__main__":
    # mpiexec -n 4 python objects/mpRfT/_2d/cf/scalar/reconstruct.py
    pass