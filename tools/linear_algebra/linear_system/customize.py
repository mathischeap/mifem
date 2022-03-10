

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
        # check 1: ---------------------------------
        if i == j:
            assert pc is None, f"when pc is None, we must have i==j, " \
                               f"now i={i}, j={j}."

        # check 2: ---------------------------------
        if pc is None:
            assert i == j, \
                f"when do not provide pc, we must set diagonal block, " \
                f"so i == j, now i={i}, j={j}."
            assert pd.__class__.__name__ == 'PartialCochain', \
                "I need a PartialCochain when pc is None."
            pc = pd

        # check 3: ---------------------------------
        assert pc.__class__.__name__ == 'PartialCochain', f"pc must be a PartialCochain."

        #======== customize ==============================================
        I, J = self._LS_.block_shape
        assert i % 1 == 0, f"i={i}({i.__class__.__name__}) cannot be an index!"
        assert j % 1 == 0, f"j={j}({j.__class__.__name__}) cannot be an index!"
        assert 0 <= i < I and 0 <= j < J, f"(i,j)= ({i},{j}) is out of range!"

        if i == j:
            self._LS_.A.customize.\
                identify_global_rows_according_to_CSCG_partial_dofs(
                i, pd, interpreted_as=interpreted_as)
            self._LS_.b.customize.\
                set_entries_according_to_two_CSCG_partial_cochains(
                i, pd, interpreted_as=interpreted_as)
        else:
            self._LS_.A.customize.\
                off_diagonally_identify_rows_according_to_two_CSCG_partial_dofs(
                i, j, pd, pc, interpreted_as=interpreted_as)
            self._LS_.b.customize.\
                set_entries_according_to_two_CSCG_partial_cochains(
                i, pd, pc=pc, interpreted_as=interpreted_as)