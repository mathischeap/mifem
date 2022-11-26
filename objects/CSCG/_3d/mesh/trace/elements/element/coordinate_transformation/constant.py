# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/08/26 2:20 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from components.freeze.base import FrozenOnly


class _3dTraceElement_CT_Constant(FrozenOnly):
    """A trace-element may have some constant CT properties. When a particular constant property
    does not exist, we return None."""

    def __init__(self, te):
        """"""
        self._te_ = te
        self._constant_unit_normal_vector_ = True # do not use None
        self._Jacobian_ = True # do not use None
        self._freeze_self_()


    @property
    def unit_normal_vector(self):
        """The direction is the right-hand-rule, it can point the inner direction of a mesh.

        When there is no constant unit_normal_vector, return None.

        Returns
        -------

        """

        if self._constant_unit_normal_vector_ is True: # not initialized as None

            te = self._te_

            if te.IS.orthogonal:
                # it has a constant_unit_normal_vector
                x = [0,]
                y = [0,]
                z = [0,]

                nV = te.coordinate_transformation.unit_normal_vector(x, y, z, parse_3_1d_eps=True)

                x, y, z = nV

                x = x[0,0]
                y = y[0,0]
                z = z[0,0]

                self._constant_unit_normal_vector_ = (x, y, z)

            # + elif, because when it is not orthogonal, the unit_normal_vector can also be constant.

            else:
                self._constant_unit_normal_vector_ = None

        return self._constant_unit_normal_vector_

    @property
    def Jacobian(self):
        """If there is no constant Jacobian for this trace element, we return None.

        Returns
        -------

        """
        if self._Jacobian_ is True: # not initialized as None

            te = self._te_

            if te.IS.orthogonal:

                tMark = te.type_wrt_metric.mark

                f0, f1 = tMark.split('d')[1:]

                f0 = float(f0[1:])
                f1 = float(f1[1:])

                self._Jacobian_ = f0 * f1 / 4


            # + elif, because when it is not orthogonal, the unit_normal_vector can also be constant.

            else:
                self._Jacobian_ = None

        return self._Jacobian_



if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
