# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/21 3:41 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
import  numpy as np


class miUsTriangular_S0F_Discretize_StandardScalar(FrozenOnly):
    """"""
    def __init__(self, sf):
        """"""
        self._sf_ = sf
        self._freeze_self_()

    def __call__(self, update_cochain=True, target='CF'):
        """
        The return cochain is 'locally full local cochain', which means it is mesh-element-wise
        local cochain. So:

        cochainLocal is a dict, whose keys are mesh element numbers, and values (1-d arrays) are
        the local cochains.
        """
        nodes = self._sf_.space.nodes
        p = self._sf_.space.p

        nodes = np.meshgrid(nodes, nodes, indexing='ij')
        nodes = [nodes[i].ravel('F')[p:] for i in range(2)]

        cochainLocal = dict()

        if target == 'CF':
            FUNC = self._sf_.CF.___DO_evaluate_func_at_time___()[0]
        else:
            raise NotImplementedError(
                f"miUsTriangular_S0F discretization does not work for target={target}.")

        for i in self._sf_.mesh.elements:
            element = self._sf_.mesh.elements[i]
            xy = element.coordinate_transformation.mapping(*nodes)
            cochainLocal[i] = FUNC(*xy)
        # pass to cochain.local ...
        if update_cochain: self._sf_.cochain.local = cochainLocal
        # ...
        return 'locally full local cochain', cochainLocal





if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
