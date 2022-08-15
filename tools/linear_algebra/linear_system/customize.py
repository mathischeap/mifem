# -*- coding: utf-8 -*-
from screws.freeze.main import FrozenOnly

class ___LinearSystem_Customize___(FrozenOnly):
    """Used to define customizations to A and b simultaneously."""
    def __init__(self, ls):
        self._LS_ = ls
        self._freeze_self_()


    def apply_strong_BC(self, i, j, pd, pc=None, interpreted_as='local_dofs'):
        """

        :param i: We apply it to block [i][j].
        :param j: We apply it to block [i][j].
        :param pd: partial degrees of freedom
        :param pc: partial cochain of freedom
        :param interpreted_as: how we interpret the `pd` and `pc`.
        :return:
        """
        if hasattr(pd, 'standard_properties') and \
             'CSCG_form' in pd.standard_properties.tags:
            pd = pd.BC.partial_cochain

        #------ for cscg meshes -------------------------------------------------------------------
        if hasattr(pd, '_mesh_') and \
            pd._mesh_.__class__.__name__ in ('_3dCSCG_Mesh', '_2dCSCG_Mesh'):

            # check 1 ____________________________________________________________________
            if i == j:
                assert pc is None, f"when i == j is None, we must have pc is None."

            # check 2 ____________________________________________________________________
            if pc is None:
                assert i == j, \
                    f"when do not provide pc, we must set diagonal block, " \
                    f"so i == j, now i={i}, j={j}."
                assert pd.__class__.__name__ == 'PartialCochain', \
                    "I need a PartialCochain when pc is None."
                pc = pd
            elif hasattr(pc, 'standard_properties') and \
                 'CSCG_form' in pc.standard_properties.tags:
                pc = pc.BC.partial_cochain
            elif hasattr(pc, 'standard_properties') and \
                 'mpRfT_form' in pc.standard_properties.tags:
                pc = pc.BC.partial_cochain
            else:
                pass

            # check 3 _______________________________________________________________________
            assert pc.__class__.__name__ == 'PartialCochain', f"pc must be a PartialCochain."

            #======== customize =============================================================
            I, J = self._LS_.block_shape
            assert i % 1 == 0, f"i={i}({i.__class__.__name__}) cannot be an index!"
            assert j % 1 == 0, f"j={j}({j.__class__.__name__}) cannot be an index!"
            assert 0 <= i < I and 0 <= j < J, f"(i,j)= ({i},{j}) is out of range!"

            if pd._mesh_.__class__.__name__ in ('_3dCSCG_Mesh', '_2dCSCG_Mesh'):
                if i == j:
                    self._LS_.A.customize.\
                        identify_global_rows_according_to_CSCG_partial_dofs(
                        i, pd, interpreted_as=interpreted_as)
                    self._LS_.b.customize.\
                        set_entries_according_to_CSCG_partial_cochains(
                        i, pd, interpreted_as=interpreted_as)
                else:
                    self._LS_.A.customize.\
                        off_diagonally_identify_rows_according_to_two_CSCG_partial_dofs(
                        i, j, pd, pc, interpreted_as=interpreted_as)
                    self._LS_.b.customize.\
                        set_entries_according_to_CSCG_partial_cochains(
                        i, pd, pc=pc, interpreted_as=interpreted_as)
            else:
                raise NotImplementedError(f"Not implemented")

        #--------- other meshes --------------------------------------------------------------------
        else:
            raise NotImplementedError(f"{pd}")


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

