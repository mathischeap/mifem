# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/10/08 10:03 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly


class miUsTriangle_SF_BC_Interpret_Local(FrozenOnly):
    """"""

    def __init__(self, sf):
        """"""
        self._sf_ = sf
        self._mesh_ = sf.mesh
        self._cochains_ = None
        self.___Pr_parse_dofs___() # need no BC.CF
        if self._sf_.BC.CF is not None:
            self.___Pr_parse_cochains___()
        else:
            pass

        self._freeze_self_()

    def ___Pr_parse_dofs___(self):
        """Need no `BC.CF`"""
        iEE = self._sf_.BC._involved_element_parts_
        self._dofs_ = dict()
        for e_e in iEE:
            element, edge = e_e
            dofs = self._sf_.numbering.do.find.local_dofs_on_element_edge(edge)
            if element not in self._dofs_:
                self._dofs_[element] = list()
            self._dofs_[element].extend(dofs)


    def ___Pr_parse_cochains___(self):
        """Need `BC.CF`"""
        indicator, local_cochain = self._sf_.discretize(target='BC')
        self._cochains_ = dict()
        if indicator == 'locally full local cochain':
            dofs = self.dofs
            for e in dofs:
                dof = dofs[e]
                cochain = local_cochain[e]
                if e not in self._cochains_:
                    self._cochains_[e] = list()

                self._cochains_[e].extend(cochain[dof])
        else:
            raise NotImplementedError()

    @property
    def dofs(self):
        return self._dofs_

    @property
    def cochains(self):
        return self._cochains_




if __name__ == "__main__":
    # mpiexec -n 4 python objects/miUsGrid/triangular/forms/standard/base/BC/interpret/local.py
    u_norm_boundaries = ['Upper', 'Down', 'Left', 'Right']
    from __init__ import miTri

    fc = miTri.call('rand0', 3)

    u = fc('1-f-o')

    u.BC.boundaries = u_norm_boundaries

    dofs = u.BC.interpret.local.dofs

    print(dofs)