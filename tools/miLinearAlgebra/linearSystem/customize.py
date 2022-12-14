# -*- coding: utf-8 -*-
from components.freeze.main import FrozenOnly

class ___LinearSystem_Customize___(FrozenOnly):
    """Used to define customizations to A and b simultaneously."""
    def __init__(self, ls):
        self._LS_ = ls
        self._freeze_self_()

    def clear(self):
        """"""
        self._LS_.A.customize.clear()
        self._LS_.b.customize.clear()

    def apply_strong_BC(self, i, j, dof_itp, cochain_itp=None, AS='local'):
        """

        :param i: We apply it to block [i][j].
        :param j: We apply it to block [i][j].
        :param dof_itp: partial degrees of freedom
        :param cochain_itp: partial cochain of freedom
        :param AS: how we interpret the boundary condition. Current it can only be
            'local`. So we apply it to local systems. It will be correctly imposed after
            assembling.
        :return:
        """
        if hasattr(dof_itp, 'standard_properties') and 'form' in dof_itp.standard_properties.tags:
            dof_itp = dof_itp.BC.interpret
        else:
            pass

        # check 1 ____________________________________________________________________
        if i == j:
            assert cochain_itp is None, f"when i == j is None, we must have pc is None."

        # check 2 ____________________________________________________________________
        if cochain_itp is None:
            assert i == j, \
                f"when do not provide pc, we must set diagonal block, " \
                f"so i == j, now i={i}, j={j}."
            cochain_itp = dof_itp

        elif hasattr(cochain_itp, 'standard_properties') and \
                'form' in cochain_itp.standard_properties.tags:
            cochain_itp = cochain_itp.BC.interpret

        else:
            pass

        #======== customize =============================================================
        I, J = self._LS_.block_shape
        assert i % 1 == 0, f"i={i}({i.__class__.__name__}) cannot be an index!"
        assert j % 1 == 0, f"j={j}({j.__class__.__name__}) cannot be an index!"
        assert 0 <= i < I and 0 <= j < J, f"(i,j)= ({i},{j}) is out of range!"
        if not isinstance(i, int): i = int(i)
        if not isinstance(j, int): j = int(j)

        if i == j:
            self._LS_.A.customize.\
                identify_global_rows_according_to(
                i, dof_itp, AS=AS)
            self._LS_.b.customize.\
                set_entries_according_to(
                i, dof_itp, AS=AS)
        else:
            self._LS_.A.customize.\
                off_diagonally_identify_rows_according_to(
                i, j, dof_itp, cochain_itp, AS=AS)
            self._LS_.b.customize.\
                set_entries_according_to(
                i, dof_itp, cochain_itp=cochain_itp, AS=AS)

    def identify_global_row(self, r):
        """We set the row #r to be all zero except M(r, r) = 1."""
        self._LS_.A.customize.identify_global_row(r)

    def set_unknown_to(self, r, v):
        """Consider ths system is Ax=b. We first clear A[r,:], then set A[r, r] = 1 and b[r] = v,
        So the unknown x[r] will equal to v.

        :param r:
        :param v:
        :return:
        """
        self._LS_.A.customize.identify_global_row(r)
        self._LS_.b.customize.set_assembled_V_i_to(r, v)