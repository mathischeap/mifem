# -*- coding: utf-8 -*-

# from components.decorators.accepts import accepts
from components.freeze.main import FrozenOnly
from scipy import sparse as spspa
import numpy as np
from root.config.main import COMM, RANK, MASTER_RANK, SIZE



class SpaVec_Customize(FrozenOnly):
    def __init__(self, spa_vec):
        self._spa_vec_ = spa_vec
        self.___customizations___ = dict()
        self._customizations_to_be_applied_ = dict()
        self._t2f_s2f_cc_ = None
        self._t1f_s1f_cc_ = None
        self._t0f_s0f_cc_ = None
        self._freeze_self_()

    def clear(self):
        """clear all existing customizations."""
        self.___customizations___ = dict()
        self._customizations_to_be_applied_ = dict()

    @property
    def _customizations_(self):
        """For example
        self.___customizations___ = {
                ...,
                5: [('cletr', 16), ...],  # The first customization for #5 element output is "clear local entry #16".
                    Customization are executed in a positive sequence.
                ...,
        }
        """
        if self.___customizations___ != dict():
            self.___PRIVATE_sent_to_applying___()
        return self._customizations_to_be_applied_


    def ___PRIVATE_sent_to_applying___(self):
        """We renew _customizations_to_be_applied_ according to ___customizations___."""

        for e in self.___customizations___:

            if e not in self._customizations_to_be_applied_:
                self._customizations_to_be_applied_[e] = (list(),  # indices
                                                          list()) # values)

            CUSi = self.___customizations___[e]  # customizations for the output of #e element.
            indices, values = self._customizations_to_be_applied_[e]

            for cus in CUSi:

                key, factors = cus


                # clear a local row of the EWC-sparse-matrix ----------------- BELOW ---------------
                if key == 'cletr':  # Clear Local EnTRy #factors
                    # `factors` will be an int
                    if factors.__class__.__name__ in ("int", "int32", "int64"):
                        pass
                    elif isinstance(factors, float):
                        assert factors % 1 == 0, f"factors={factors} wrong, when `cletr`, factors must be an int."
                    else:
                        raise Exception(
                            f"factors={factors} wrong, when `cletr`, factors must be an int.")

                    if factors not in indices:
                        indices.append(factors)
                        values.append(0)
                    else:
                        ind = indices.index(factors)
                        values[ind] = 0


                elif key == 'slet':  # set local entry to
                    i, v = factors

                    if i not in indices:
                        indices.append(i)
                        values.append(v)
                    else:
                        ind = indices.index(i)
                        values[ind] = v


                elif key == 'slest':  # Set Local Entries To
                    I, V = factors
                    for i, v in zip(I, V):
                        if i not in indices:
                            indices.append(i)
                            values.append(v)
                        else:
                            ind = indices.index(i)
                            values[ind] = v

                # Not Implemented ---------------- BELOW ----------------------------
                else:
                    raise NotImplementedError(f"Can not handle customization key={key}.")
                # ================================ ABOVE ============================

        self.___customizations___ = dict() # we have to clear ___customizations___ to avoid multiple renewing.


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
            RETURN = RETURN.tolil()
            indices, values = self._customizations_[e]
            RETURN[indices, 0] = values

        if not spspa.isspmatrix_csc(RETURN): RETURN = RETURN.tocsc()
        assert RETURN.shape[1] == 1, f"A vector must to csc_matrix of shape (x, 1), " \
                                     f"now its shape is {RETURN.shape}."
        return RETURN


    def set_entries_according_to(self, i, BC_itp, cochain_itp=None, AS='local'):
        """

        Parameters
        ----------
        i
        BC_itp
        cochain_itp
        AS

        Returns
        -------

        """
        assert self._spa_vec_.elements._mesh_ == BC_itp._mesh_, "BC_itp mesh != EWC mesh."
        if cochain_itp is not None:
            assert self._spa_vec_.elements._mesh_ == cochain_itp._mesh_, "cochain_itp mesh != EWC mesh."

        bsp = self._spa_vec_.con_shape
        assert bsp is not False and np.shape(bsp) == (1,), \
            "we can only use this to concatenate EWC vectors to keep the block safe. " \
            "if you want to apply it to a single EWC vectors, first wrap it in a list, then send it to concatenate"

        assert i % 1 == 0, f"i={i}({i.__class__.__name__}) cannot be an index!"
        if isinstance(i, float): i = int(i)
        assert 0 <= i < bsp[0], f"i={i} is beyond the shape of this " \
                                f"EWC spare vector of block shape {bsp}."


        if AS == 'local':
            dofs = BC_itp.local.dofs
            if cochain_itp is None:
                cochains = BC_itp.local.cochains
            else:
                cochains = cochain_itp.local.cochains


            GM = self._spa_vec_.gathering_matrix
            if GM.___Pr_IS_regular___:
                LDR_row = GM.___Pr_regular__local_dofs_ranges___[i]
                start_row = LDR_row.start
            else:
                raise NotImplementedError(f"Not implemented for irregular Gathering Matrix.")

            for e in dofs:

                local_dofs = np.array(dofs[e], dtype=int) + start_row
                ENTRIES = cochains[e]

                assert local_dofs.shape[0] == len(ENTRIES)

                if e not in self.___customizations___:
                    self.___customizations___[e] = list()

                self.___customizations___[e].append(('slest', (local_dofs, ENTRIES)))

        else:
            raise NotImplementedError()

    def set_constant_entries_according_to(self, i, BC_itp, constant, AS='local'):
        """

        Parameters
        ----------
        i
        BC_itp :
            For example, `u.BC.interpret`.
        constant
        AS

        Returns
        -------

        """
        assert self._spa_vec_.elements._mesh_ == BC_itp._mesh_, "BC_itp mesh != EWC mesh."

        bsp = self._spa_vec_.con_shape
        assert bsp is not False and np.shape(bsp) == (1,), \
            "we can only use this to concatenate EWC vectors to keep the block safe. " \
            "if you want to apply it to a single EWC vectors, first wrap it in a list, then send it to concatenate"

        assert i % 1 == 0, f"i={i}({i.__class__.__name__}) cannot be an index!"
        if isinstance(i, float): i = int(i)
        assert 0 <= i < bsp[0], f"i={i} is beyond the shape of this " \
                                f"EWC spare vector of block shape {bsp}."


        if AS == 'local':
            dofs = BC_itp.local.dofs

            GM = self._spa_vec_.gathering_matrix
            if GM.___Pr_IS_regular___:
                LDR_row = GM.___Pr_regular__local_dofs_ranges___[i]
                start_row = LDR_row.start
            else:
                raise NotImplementedError(f"Not implemented for irregular Gathering Matrix.")

            for e in dofs:

                local_dofs = np.array(dofs[e], dtype=int) + start_row
                CONSTANT = constant * np.ones(len(local_dofs))

                if e not in self.___customizations___:
                    self.___customizations___[e] = list()

                self.___customizations___[e].append(('slest', (local_dofs, CONSTANT)))


        else:
            raise NotImplementedError(f"interpret={AS} not implemented.")

    #
    # @accepts('self', (int, float, 'int32', 'int64'), ('PartialDofs', 'PartialCochain'))
    # def set_entries_according_to_CSCG_partial_cochains(
    #     self, i, pd, pc=None, interpreted_as='local_dofs'):
    #     """We do:
    #     block[i][row_pds.dofs, 1] = col_pds.cochain
    #     locally. Note it is "locally", that means, after assembly, we
    #     may have values which are multiple times of the value. This is okay because same thing
    #     will happen in the left-hand-side EWC matrix.
    #
    #     :param i:
    #     :param pd:
    #     :param pc:
    #     :param interpreted_as:
    #     :return:
    #     """
    #     if pc is None:
    #         assert pd.__class__.__name__ == 'PartialCochain', \
    #             "I need a PartialCochain when pc is None."
    #         pc = pd
    #
    #     if pd.__class__.__name__ == 'PartialCochain': pd = pd.dofs
    #
    #     assert pd.__class__.__name__ == 'PartialDofs', f"pd must be a PartialDofs."
    #     assert pc.__class__.__name__ == 'PartialCochain', f"pc must be a PartialCochain."
    #
    #     bsp = self._spa_vec_.con_shape
    #     assert bsp is not False and np.shape(bsp) == (1,), \
    #         "we can only use this to concatenate EWC vectors to keep the block safe. " \
    #         "if you want to apply it to a single EWC vectors, first wrap it in a list, then send it to concatenate"
    #
    #     assert i % 1 == 0, f"i={i}({i.__class__.__name__}) cannot be an index"
    #     if isinstance(i, float): i = int(i)
    #     assert 0 <= i < bsp[0], f"i={i} is beyond the shape of this " \
    #                             f"EWC spare vector of block shape {bsp}."
    #
    #
    #     if interpreted_as == 'local_dofs':
    #         GM = self._spa_vec_.gathering_matrix
    #         if GM.___Pr_IS_regular___:
    #             LDR_row = GM.___Pr_regular__local_dofs_ranges___[i]
    #             start_row = LDR_row.start
    #         else:
    #             raise NotImplementedError(f"Not implemented for irregular Gathering Matrix.")
    #         LDF_row = pd.interpreted_as.local_dofs
    #         # LDF_col = pc.dofs.interpreted_as.local_dofs
    #         cochain = pc.cochain
    #         assert len(LDF_row) == len(cochain), "pc, pd must cover same amount of mesh elements."
    #         for e in LDF_row:
    #
    #             ROW = LDF_row[e]
    #             # COL = LDF_col[e]
    #             LOC = cochain[e]
    #
    #             # LOC = self.___PRIVATE_correcting_correspondence___(pd._form_, ROW, COL, LOC, pc._form_)
    #
    #             local_dofs = np.array(ROW) + start_row
    #             local_cochain = LOC
    #
    #             if e not in self.___customizations___:
    #                 self.___customizations___[e] = list()
    #             # Set Local Entries
    #             self.___customizations___[e].append(('slest', (local_dofs, local_cochain)))
    #
    #     else:
    #         raise NotImplementedError(f"interpreted_as={interpreted_as} not implemented.")
    #
    #
    # # def ___PRIVATE_correcting_correspondence___(self, rf, R, C, local_cochain, cf):
    # #     """
    # #
    # #     Parameters
    # #     ----------
    # #     rf
    # #     R
    # #     C
    # #     cf
    # #
    # #     Returns
    # #     -------
    # #
    # #     """
    # #     if rf.__class__.__name__ == '_3dCSCG_1Trace' and cf.__class__.__name__ == '_3dCSCG_1Form':
    # #         # we have to make sure this to make the singularity handling possible
    # #         if self._t1f_s1f_cc_ is None:
    # #             s1f = list()
    # #             for side in 'NSWEBF':
    # #                 s1f.append(cf.numbering.do.find.local_dofs_on_element_side(side))
    # #             s1f = np.concatenate(s1f)
    # #
    # #             self._t1f_s1f_cc_ = s1f
    # #         else:
    # #             s1f = self._t1f_s1f_cc_
    # #
    # #         cC = s1f[R]
    # #         assert len(cC) == len(C) and set(cC) == set(C), f"must be this case."
    # #
    # #         COCHAIN = np.empty(cf.num.basis)
    # #         COCHAIN[C] = local_cochain
    # #         return COCHAIN[cC]
    # #
    # #     elif rf.__class__.__name__ == '_3dCSCG_0Trace' and cf.__class__.__name__ == '_3dCSCG_0Form':
    # #         # we have to make sure this to make the singularity handling possible
    # #         if self._t0f_s0f_cc_ is None:
    # #             s0f = list()
    # #             for side in 'NSWEBF':
    # #                 s0f.append(cf.numbering.do.find.local_dofs_on_element_side(side))
    # #             s0f = np.concatenate(s0f)
    # #
    # #             self._t0f_s0f_cc_ = s0f
    # #         else:
    # #             s0f = self._t0f_s0f_cc_
    # #
    # #         cC = s0f[R]
    # #         assert len(cC) == len(C) and set(cC) == set(C), f"must be this case."
    # #
    # #         COCHAIN = np.empty(cf.num.basis)
    # #         COCHAIN[C] = local_cochain
    # #         return COCHAIN[cC]
    # #
    # #     else:
    # #         return local_cochain
    #
    #
    #










    def set_assembled_V_i_to(self, i, v):
        """Let V be the assembled vector (csc_matrix of shape (x,1)), we set V[i] = v.

        :param i:
        :param v:
        :return:
        """
        assert i % 1 == 0, f"i={i}({i.__class__.__name__}) is wrong."
        assert isinstance(v, (int, float))

        if not isinstance(i, int): i = int(i)

        rCGM = self._spa_vec_.gathering_matrix
        assert rCGM is not None, "I have no gathering_matrix!"

        rO = rCGM.do.find.elements_and_local_indices_of_dof(i)

        if rO is None:
            elements, indices = None, None
        else:
            elements, indices = rO

        all_elements = COMM.gather(elements, root=MASTER_RANK)

        if RANK == MASTER_RANK:
            ELEMENTS = set()
            for _ES_ in all_elements:
                if _ES_ is not None:
                    ELEMENTS.update(_ES_)

            ELEMENTS = list(ELEMENTS)
            assert len(ELEMENTS) > 0, f"I do not find any!"
            ELEMENTS.sort()
            the_element = ELEMENTS[0]
        else:
            the_element = None
        the_element = COMM.bcast(the_element, root=MASTER_RANK)

        if rO is None:
            I_am_in = SIZE
        else:
            if the_element in elements:
                I_am_in = RANK
            else:
                I_am_in = SIZE

        who_are_in = COMM.gather(I_am_in, root=MASTER_RANK)
        if RANK == MASTER_RANK:
            use_who = min(who_are_in)
            assert use_who != SIZE, "Something is wrong, we need to find a correct core."
        else:
            use_who = None
        use_who = COMM.bcast(use_who, root=MASTER_RANK)

        if use_who == RANK:
            assert the_element == elements[0], "Must be like this!"

            for e, i in zip(elements, indices):

                if e not in self.___customizations___:
                    self.___customizations___[e] = list()

                if e == the_element:
                    self.___customizations___[e].append(('slet', (i, v)))
                else:
                    self.___customizations___[e].append(('cletr', i))

        else:
            if elements is None:
                pass
            else:
                for e, i in zip(elements, indices):
                    if e not in self.___customizations___:
                        self.___customizations___[e] = list()
                    self.___customizations___[e].append(('cletr', i))