# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/10/13 9:19 PM
"""

from components.freeze.base import FrozenOnly


class CSCG_AFORM_BC_Interpret_Local(FrozenOnly):
    """"""

    def __init__(self, adf):
        """"""
        self._adf_ = adf
        self._mesh_ = adf.mesh
        self._cochains_ = None
        self.___Pr_parse_dofs___()  # need no BC.CF
        if self._adf_.BC.CF is not None:
            self.___Pr_parse_cochains___()
        else:
            pass
        self._freeze_self_()

    def ___Pr_parse_dofs___(self):
        """Need no `BC.CF`"""
        iEP = self._adf_.BC._involved_element_parts_
        self._dofs_ = dict()
        for e_p in iEP:
            element, part = int(e_p[:-1]), e_p[-1]
            if self._mesh_.ndim == 3:
                dofs = self._adf_.prime.numbering.do.\
                    find.local_dofs_on_element_side(part)
            elif self._mesh_.ndim == 2:
                dofs = self._adf_.prime.numbering.do.\
                    find.local_dofs_on_element_edge(part)
            else:
                raise Exception()
            if element not in self._dofs_:
                self._dofs_[element] = list()
            self._dofs_[element].extend(dofs)

    def ___Pr_parse_cochains___(self):
        """Need `BC.CF`"""
        indicator, local_cochain = self._adf_.prime.discretize(target='BC')
        self._cochains_ = dict()

        if indicator == 'locally full local cochain':

            raise NotImplementedError()

            # dofs = self.dofs
            # for e in dofs:
            #     dof = dofs[e]
            #     cochain = local_cochain[e]
            #     if e not in self._cochains_:
            #         self._cochains_[e] = list()
            #
            #     self._cochains_[e].extend(cochain[dof])

        elif indicator == 'locally full local TEW cochain':
            
            iEP = self._adf_.BC._involved_element_parts_
            t = self._adf_  # must be an adf trace form.
            TEM = t.mesh.trace.elements.map

            MM_TEW = t.prime.___PRIVATE_generate_TEW_mass_matrices___()

            if t.ndim == 2:
                side_index = {'U': 0, 'D': 1, 'L': 2, 'R': 3}
            elif t.ndim == 3:
                side_index = {'N': 0, 'S': 1, 'W': 2, 'E': 3, 'B': 4, 'F': 5}
            else:
                raise Exception()

            for e_p in iEP:
                element, part = int(e_p[:-1]), e_p[-1]

                if element not in self._cochains_:
                    self._cochains_[element] = list()
                index = side_index[part]
                TE = TEM[element][index]
                prime_cc_TEW = local_cochain[TE]
                dual_cc_TEW = MM_TEW[TE] @ prime_cc_TEW
                self._cochains_[element].extend(dual_cc_TEW)

        else:
            raise NotImplementedError(f"{indicator}")

    @property
    def dofs(self):
        return self._dofs_

    @property
    def cochains(self):
        return self._cochains_
