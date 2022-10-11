# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/21 4:19 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from objects.miUsGrid.base.form.standard.cochain import miUsGrid_SF_CochainBase

class miUs_Triangular_SF_Cochain(miUsGrid_SF_CochainBase):
    """"""

    def __init__(self, sf):
        """"""
        super(miUs_Triangular_SF_Cochain, self).__init__(sf)
        self._freeze_self_()

    def ___PRIVATE_local_on_axis___(self, axis, i):
        """
        find the local cochain along a particular axis in a particular mesh_element

        :param str axis: The local cochain along which axis? ``x``, ``y`` or ``z``.
        :param int i: in this mesh element
        :return: The local cochain dict.
        :rtype: Dict[int, numpy.ndarray]
        """
        numOfBasisComponents = getattr(self._sf_.space.num_basis_components,
                                       self._sf_.__class__.__name__)
        if axis == 'x':
            localAlongAxis = self.local[i][:numOfBasisComponents[0]]
        elif axis == 'y':
            localAlongAxis = self.local[i][numOfBasisComponents[0]:]
        else:
            raise Exception()
        return localAlongAxis



if __name__ == "__main__":
    # mpiexec -n 4 python objects/miUsGrid/triangular/forms/standard/base/cochain/main.py
    pass
