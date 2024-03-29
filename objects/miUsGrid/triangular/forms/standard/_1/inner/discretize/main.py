# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/21 6:24 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from components.freeze.base import FrozenOnly
from objects.miUsGrid.triangular.forms.standard._1.inner.discretize.vector.standard import \
    miUsTriangular_iS1F_Discretize_StandardVector


class miUsTriangular_iS1F_Discretize(FrozenOnly):
    def __init__(self, sf):
        self._sf_ = sf
        self._standard_vector_ = miUsTriangular_iS1F_Discretize_StandardVector(sf)
        self._freeze_self_()

    def __call__(self, target='CF', **kwargs):
        """
        Discretize the current function (a scalar field) to cochain.

        It is actually a wrapper of multiple methods that discretize functions of different types.

        """
        if target == 'CF':

            if self._sf_.CF.__class__.__name__ == 'miUsGrid_Triangular_Vector':
                #
                if self._sf_.CF.ftype == 'standard':
                    return self._standard_vector_(**kwargs)
                else:
                    raise NotImplementedError(f"miUsTriangular_iS1F cannot (target func) "
                                              f"discretize miUsGrid_Triangular_Vector of ftype={self._sf_.CF.ftype}")

            else:
                raise NotImplementedError(f'miUsTriangular_iS1F can not (target func) '
                                          f'discretize {self._sf_.CF.__class__}.')

        elif target == 'BC':

            if self._sf_.BC.CF.__class__.__name__ == 'miUsGrid_Triangular_Vector':

                if self._sf_.BC.CF.ftype == 'standard':
                    return self._standard_vector_(update_cochain=False, target='BC', **kwargs)

                else:
                    raise NotImplementedError()

            else:
                raise NotImplementedError()


        else:
            raise NotImplementedError(f"miUsTriangular_iS1F cannot discretize "
                                      f"while targeting at {target}.")

if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
