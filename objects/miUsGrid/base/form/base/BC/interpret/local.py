# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/10/08 10:03 PM
"""
import sys
if './' not in sys.path: sys.path.append('./')
from components.freeze.base import FrozenOnly

# look, we only cache the dofs, not the cochains.
___miUSGrid_global_BC_itp_Local_dofs_cache___ = {
    'iEE_ID': -1,
    'dofs': dict(),
}

class miUsGrid_Form_BC_Interpret_Local(FrozenOnly):
    """"""

    def __init__(self, f):
        """"""
        self._f_ = f
        self._mesh_ = f.mesh
        self._cochains_ = None
        self.___Pr_parse_dofs___() # need no BC.CF
        if self._f_.BC.CF is not None:
            self.___Pr_parse_cochains___()
        else:
            pass
        self._freeze_self_()

    def ___Pr_parse_dofs___(self):
        """Need no `BC.CF`"""
        iEE = self._f_.BC._involved_element_parts_

        if id(iEE) == ___miUSGrid_global_BC_itp_Local_dofs_cache___['iEE_ID']:
            self._dofs_ = ___miUSGrid_global_BC_itp_Local_dofs_cache___['dofs']

        else:
            self._dofs_ = dict()
            for e_e in iEE:
                element, edge = e_e

                if self._mesh_.ndim == 3:
                    raise NotImplementedError()
                elif  self._mesh_.ndim == 2:
                    dofs = self._f_.numbering.do.find.local_dofs_on_element_edge(edge)
                else:
                    raise Exception()

                if element not in self._dofs_:
                    self._dofs_[element] = list()
                self._dofs_[element].extend(dofs)

            ___miUSGrid_global_BC_itp_Local_dofs_cache___['iEE_ID'] = id(iEE)
            ___miUSGrid_global_BC_itp_Local_dofs_cache___['dofs'] = self._dofs_

    def ___Pr_parse_cochains___(self):
        """Need `BC.CF`"""
        indicator, local_cochain = self._f_.discretize(target='BC')
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

    def __iter__(self):
        for e in self.dofs:
            yield e

    def __len__(self):
        return len(self.dofs)

    def __contains__(self, item):
        return item in self.dofs

    @property
    def dofs(self):
        return self._dofs_

    @property
    def cochains(self):
        return self._cochains_

if __name__ == "__main__":
    # mpiexec -n 4 python objects/miUsGrid/base/form/base/BC/interpret/local.py
    u_norm_boundaries = ['Upper', 'Down', 'Left', 'Right']
    from __init__ import miTri

    fc = miTri.call('rand0', 3)

    u = fc('1-f-o')

    u.BC.boundaries = u_norm_boundaries

    dofs = u.BC.interpret.local.dofs

    print(dofs)