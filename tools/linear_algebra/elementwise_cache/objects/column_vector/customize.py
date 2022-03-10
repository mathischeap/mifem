


from screws.decorators.accepts import accepts
from screws.freeze.main import FrozenOnly
from scipy import sparse as spspa
import numpy as np



class SpaVec_Customize(FrozenOnly):
    def __init__(self, spa_vec):
        self._spa_vec_ = spa_vec
        self.___customizations___ = dict()
        self._freeze_self_()

    @property
    def _customizations_(self):
        """For example
        self.___customizations___ = {
                ...,
                5: [('cletr', 16), ...],  # The first customization for #5 element output is "clear local entry #16".
                    Customization are executed in a positive sequence.
                ...
        }
        """
        return self.___customizations___

    def ___PRIVATE_do_execute_customization___(self, RETURN, e):
        """Execute all added customization.

        :param RETURN: The RETURN to be customized.
        :param e: for #i element.
        :return:
        """
        if e not in self._customizations_:

            pass

        else:
            assert spspa.isspmatrix_csc(RETURN), "We need to start with a (1d) csc matrix."

            CUSi = self._customizations_[e] # customizations for the output of #e element.
            for cus in CUSi:

                key, factors = cus

                # clear a local row of the EWC-sparse-matrix ----------------- BELOW -----------------------------------
                if key == 'cletr': # Clear Local EnTRy #factors
                    # `factors` will be an int
                    if factors.__class__.__name__ in ("int", "int32", "int64"):
                        pass
                    elif isinstance(factors, float):
                        assert factors % 1 == 0, f"factors={factors} wrong, when `cletr`, factors must be an int."
                    else:
                        raise Exception(f"factors={factors} wrong, when `cletr`, factors must be an int.")

                    if not spspa.isspmatrix_lil(RETURN): RETURN = RETURN.tolil()
                    RETURN[factors, 0] = 0
                elif key == 'slest': # Set Local Entries To
                    if not spspa.isspmatrix_lil( RETURN): RETURN = RETURN.tolil()
                    RETURN[factors[0], 0] = factors[1]

                # Not Implemented --------------------------------------------- BELOW ----------------------------------
                else:
                    raise NotImplementedError(f"Can not handle customization key={key}.")
                #============================================================== ABOVE ==================================

        if not spspa.isspmatrix_csc(RETURN): RETURN = RETURN.tocsc()
        assert RETURN.shape[1] == 1, f"A vector must to csc_matrix of shape (x, 1), now its shape is {RETURN.shape}."
        return RETURN

    @accepts('self', (int, float, 'int32', 'int64'),
             ('PartialDofs', 'PartialCochain'), (int, float))
    def set_constant_entries_according_to_CSCG_partial_dofs(
            self, i, pd, constant, interpreted_as='local_dofs'):
        """ We do:
        block[i][row_pds.dofs, 1] = constant
        locally. Note it is "locally", that means, after assembly, we
        may have values which are multiple times of the constant. This is okay because same thing
        will happen in the left-hand-side EWC matrix.

        :param i:
        :param pd:
        :param constant:
        :param interpreted_as:
        :return:
        """

        if pd.__class__.__name__ == 'PartialCochain': pd = pd.dofs

        bsp = self._spa_vec_.con_shape
        assert bsp is not False and np.shape(bsp) == (1,), \
            "we can only use this to concatenate EWC vectors to keep the block safe. " \
            "if you want to apply it to a single EWC vectors, first wrap it in a list, then send it to concatenate"

        assert i % 1 == 0, f"i={i}({i.__class__.__name__}) cannot be an index!"
        if isinstance(i, float): i = int(i)
        assert 0 <= i < bsp[0], f"i={i} is beyond the shape of this " \
                                f"EWC spare vector of block shape {bsp}."

        if interpreted_as == 'local_dofs':
            GM = self._spa_vec_.gathering_matrix
            LDR_row = GM.local_dofs_ranges[i]
            start_row = LDR_row.start
            dofs = pd.interpreted_as.local_dofs
            for e in dofs:
                local_dofs = np.array(dofs[e]) + start_row

                if e not in self.___customizations___:
                    self.___customizations___[e] = list()

                self.___customizations___[e].append(('slest', (local_dofs, constant)))

        else:
            raise NotImplementedError(f"interpreted_as={interpreted_as} not implemented.")

    @accepts('self', (int, float, 'int32', 'int64'), ('PartialDofs', 'PartialCochain'))
    def set_entries_according_to_two_CSCG_partial_cochains(
        self, i, pd, pc=None, interpreted_as='local_dofs'):
        """We do:
        block[i][row_pds.dofs, 1] = col_pds.cochain
        locally. Note it is "locally", that means, after assembly, we
        may have values which are multiple times of the value. This is okay because same thing
        will happen in the left-hand-side EWC matrix.

        :param i:
        :param pd:
        :param pc:
        :param interpreted_as:
        :return:
        """
        if pc is None:
            assert pd.__class__.__name__ == 'PartialCochain', \
                "I need a PartialCochain when pc is None."
            pc = pd

        if pd.__class__.__name__ == 'PartialCochain': pd = pd.dofs

        assert pd.__class__.__name__ == 'PartialDofs', f"pd must be a PartialDofs."
        assert pc.__class__.__name__ == 'PartialCochain', f"pc must be a PartialCochain."

        bsp = self._spa_vec_.con_shape
        assert bsp is not False and np.shape(bsp) == (1,), \
            "we can only use this to concatenate EWC vectors to keep the block safe. " \
            "if you want to apply it to a single EWC vectors, first wrap it in a list, then send it to concatenate"

        assert i % 1 == 0, f"i={i}({i.__class__.__name__}) cannot be an index!"
        if isinstance(i, float): i = int(i)
        assert 0 <= i < bsp[0], f"i={i} is beyond the shape of this " \
                                f"EWC spare vector of block shape {bsp}."


        if interpreted_as == 'local_dofs':
            GM = self._spa_vec_.gathering_matrix
            LDR_row = GM.local_dofs_ranges[i]
            start_row = LDR_row.start
            dofs = pd.interpreted_as.local_dofs
            cochain = pc.cochain
            assert len(dofs) == len(cochain), "pc, pd must cover same amount of mesh elements."
            for e in dofs:
                local_dofs = np.array(dofs[e]) + start_row
                local_cochain = cochain[e]
                if e not in self.___customizations___:
                    self.___customizations___[e] = list()
                # Set Local Entries To
                self.___customizations___[e].append(('slest', (local_dofs, local_cochain)))

        else:
            raise NotImplementedError(f"interpreted_as={interpreted_as} not implemented.")