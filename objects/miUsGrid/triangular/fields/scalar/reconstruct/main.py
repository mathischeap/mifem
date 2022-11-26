# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/21 2:14 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from components.freeze.base import FrozenOnly
from objects.miUsGrid.triangular.fields.scalar.reconstruct.mesh_element.standard import OnMeshElement_for_Standard

class miUsGrid_Triangular_Scalar_Reconstruct(FrozenOnly):
    """"""

    def __init__(self, scalar):
        """"""
        self._scalar_ = scalar
        self._on_mesh_element___for_standard_ = OnMeshElement_for_Standard(scalar)
        self._freeze_self_()

    def __call__(self, xi, eta, time=None, ravel=False, i=None, where=None):
        """

        :param xi: We will map (xi, eta) to physical domain without any processing like `meshgrid`.
        :param eta: We will map (xi, eta) to physical domain without any processing like `meshgrid`.
        :param time:
        :param ravel:
        :param i:
        :param where:
        :return:
        """
        ftype = self._scalar_.ftype

        #------- deal with time --------------------------------------------------
        if time is None:
            pass
        else:
            self._scalar_.current_time = time

        # parse where when it is None --------------------------------------------
        if where is None:
            if ftype == 'standard':
                where = 'mesh-element'
            else:
                raise NotImplementedError(f"please set default `where` for {ftype} TriScalarField.")
        else:
            pass

        #--------------------------------------------------------------------------
        if where == 'mesh-element':
            if ftype == 'standard':
                return self._on_mesh_element___for_standard_(xi, eta, ravel, i)
            else:
                raise NotImplementedError()
        else:
            raise NotImplementedError()


if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
