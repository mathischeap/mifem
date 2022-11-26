# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/10/13 9:19 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from components.freeze.base import FrozenOnly


class CSCG_FORM_BC_Interpret_Local(FrozenOnly):
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
        iEP = self._f_.BC._involved_element_parts_
        self._dofs_ = dict()
        for e_p in iEP:
            element, part = int(e_p[:-1]), e_p[-1]
            if  self._mesh_.ndim == 3:
                dofs = self._f_.numbering.do.\
                    find.local_dofs_on_element_side(part)
            elif self._mesh_.ndim == 2:
                dofs = self._f_.numbering.do.\
                    find.local_dofs_on_element_edge(part)
            else:
                raise Exception()
            if element not in self._dofs_:
                self._dofs_[element] = list()
            self._dofs_[element].extend(dofs)

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

        elif indicator == 'Boundary only local cochain':
            # only cochains for dofs on mesh element side (boundary of the mesh).
            iEP = self._f_.BC._involved_element_parts_

            for e_p in iEP:
                e, part = int(e_p[:-1]), e_p[-1]
                assert e in local_cochain, \
                    f"element {e} in not in the local cochain, most likely," \
                    f"the boundaries in the func do not cover BC.valid_boundaries: {self._f_.BC.boundaries}."

                if e not in self._cochains_:
                    self._cochains_[e] = list()

                assert len(self.dofs[e]) > 0, f"empty for element #{e}"

                self._cochains_[e].extend(local_cochain[e][part])


        elif indicator == 'locally full local TEW cochain':

            iEP = self._f_.BC._involved_element_parts_
            t = self._f_ # must be a trace form.
            TEM = t.mesh.trace.elements.map
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
                self._cochains_[element].extend(local_cochain[TE])

        else:
            raise NotImplementedError(f"{indicator}")

    @property
    def dofs(self):
        return self._dofs_

    @property
    def cochains(self):
        return self._cochains_





if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
