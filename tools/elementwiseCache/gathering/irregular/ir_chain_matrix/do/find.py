# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/06/12 5:02 PM
"""
import sys

if './' not in sys.path: sys.path.append('/')

from components.freeze.base import FrozenOnly
from functools import lru_cache
from components.decorators.all import accepts
import numpy as np


class iR_CGM_DO_FIND(FrozenOnly):
    """"""

    def __init__(self, ir_cgm):
        """"""
        self._CGM_ = ir_cgm
        self._freeze_self_()

    @accepts('self', (int, float, 'int32', 'int64'))
    def elements_contain_dof_numbered(self, m, N=None):
        """Find the local mesh-element(s) which contains the dof numbered `m`.

        :param m: the dof to be found.
        :param N: The max number of elements to find. If we already find `N` elements containing the dof `m`, we stop
            finding. When `N=None`, then we always check all local elements.
        :return:

            A list of str(s) represent the local element(s) that has the dof numbered `m` or None if no
            local element contains that dof. Therefore, maybe in multiple cores, the method does not return None.

            We also return `where` (int): GMs[where] has the dof `m`. When the first output is None, we do not return where.
        """
        assert (m % 1) == 0, f"m={m} is wrong."
        assert -self._CGM_.global_num_dofs <= m < self._CGM_.global_num_dofs, \
                f"dof numbered = {m} is out of range, it should be in " \
                f"[{-self._CGM_.global_num_dofs}, {self._CGM_.global_num_dofs - 1}] " \
                f"cause the maximum global numbering for this chain_gathering_matrix " \
                f"is {self._CGM_.global_num_dofs - 1}."

        if m < 0: m += self._CGM_.global_num_dofs
        assert 0 <= m < self._CGM_.global_num_dofs, f"m={m} is wrong."
        if not isinstance(m, int): m = int(m)

        if N is None:
            pass
        else:
            pass

        assert N is None or N > 0, f"At least we search for 1 element, right? Now it is N={N}"

        if len(self._CGM_) == 0: # iR_chain_gathering_matrix is empty: no local element at all.
            return None

        if self._CGM_.chain_method == 'silly':
            LRS = self._CGM_.local_ranges

            where = -1

            if self._CGM_.num_GMs == 1:
                MIN, MAX = LRS[0]
                if MIN <= m < MAX:
                    where = 0
                else:
                    return None

            else:
                for i, LR in enumerate(LRS):
                    MIN, MAX = LR
                    if i == 0: # the first GM
                        if m < MIN:
                            return None
                        elif MIN <= m < MAX:
                            where = i
                            break
                        else:
                            pass

                    else:

                        if LRS[i-1][1] <= m < MIN:
                            return None
                        elif MIN <= m < MAX:
                            where = i
                            break
                        else:
                            if i == self._CGM_.num_GMs - 1: # the last GM
                                if m >= MAX:
                                    return None
                            else:
                                pass

            # if we reach here, we know this core may have element(s) containing the target dof...
            assert 0 <= where < self._CGM_.num_GMs, "MUST BE!" # we know we only need to look at GMs[where] to find m.

            m -= self._CGM_._To_Be_Added_[where]
            GM = self._CGM_.GMs[where]

            # now we go through `GM` to find the element(s) having `m`
            ELE = list()

            if N is None:
                for rp in GM:
                    if m in GM[rp]: # NICE! we find one element containing m.
                        ELE.append(rp)

            else:
                n = 0
                for rp in GM:
                    if m in GM[rp]: # NICE! we find one element containing m.
                        ELE.append(rp)
                        n += 1

                        if n == N: # we have found enough elements. Let's break the loop.
                            break

            if ELE == list():
                return None
            else:
                return ELE, int(where)

        else:
            raise NotImplementedError()


    @accepts('self', (int, float, 'int32', 'int64'))
    @lru_cache(maxsize=256)
    def elements_and_local_indices_of_dof(self, m, N=None):
        """
        Find the local element(s) and local indices of the dof numbered i.

        :param m:
        :param N:
        :return: If no element contains dof m, return None. Else return two lists, one is for element(s), one is
            for the corresponding local index(s).
        """
        assert (m % 1) == 0, f"m={m} is wrong."
        assert -self._CGM_.global_num_dofs <= m < self._CGM_.global_num_dofs, \
                f"dof numbered = {m} is out of range, it should be in " \
                f"[{-self._CGM_.global_num_dofs}, {self._CGM_.global_num_dofs - 1}] " \
                f"cause the maximum global numbering for this chain_gathering_matrix " \
                f"is {self._CGM_.global_num_dofs - 1}."
        if m < 0: m += self._CGM_.global_num_dofs
        assert 0 <= m < self._CGM_.global_num_dofs, f"m={m} is wrong."

        if not isinstance(m, int): m = int(m)

        OUT = self.elements_contain_dof_numbered(m, N=N)

        if self._CGM_.chain_method == 'silly':
            if OUT is None:
                return None
            elements = OUT[0]
            local_indices = list()

            for rp in elements:

                index = np.argwhere(self._CGM_[rp]==m)[0,0]
                local_indices.append(int(index))

            return elements, local_indices

        else:
            mesh_elements = list()
            local_indices = list()

            for e in self._CGM_: # go through all local mesh elements
                gv = self._CGM_[e] # get the local gathering vector in each local mesh element
                if m in gv:
                    mesh_elements.append(e)
                    local_indices.append(gv.index(m))

            if len(mesh_elements) == 0:
                return None
            else:
                return mesh_elements, local_indices

    def elements_and_local_indices_of_dofs(self, dofs):
        """"""
        mesh_elements = dict()
        local_indices = dict()

        for i in dofs:
            mesh_elements[i] = list()
            local_indices[i] = list()

        for e in self._CGM_: # go through all local mesh elements
            gv = self._CGM_[e] # get the local gathering vector in each local mesh element
            for i in dofs:
                if i in gv:
                    mesh_elements[i].append(e)
                    local_indices[i].append(gv.index(i))

        return mesh_elements, local_indices



if __name__ == "__main__":
    # mpiexec -n 4 python tools/linear_algebra/gathering/irregular/ir_chain_matrix/do/find.py
    pass
