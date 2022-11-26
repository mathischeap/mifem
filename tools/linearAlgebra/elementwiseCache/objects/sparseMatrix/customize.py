# -*- coding: utf-8 -*-

from components.freeze.main import FrozenOnly
from scipy import sparse as spspa
from root.config.main import *


class SpaMat_Customize(FrozenOnly):
    """Store all method to customize the EWC sparse matrix.

    Notice that all customization will only affect at the level of this SpaMat, its blocks will
    not be affected. So if we make other SpaMat with its blocks, those blocks will remain
    not-customized.

    """
    def __init__(self, spa_mat):
        self._spa_mat_ = spa_mat
        self._t2f_s2f_cc_ = None
        self._t1f_s1f_cc_ = None
        self._t0f_s0f_cc_ = None
        self.___customizations___ = dict()
        self._freeze_self_()

    def clear(self):
        """clear all existing customizations."""
        self.___customizations___ = dict()

    @property
    def _customizations_(self):
        """For example
        self.___customizations___ = {
                ...,
                5: [('clr', 16), ...],  # The first customization for #5 element output is "clear local row #16".
                    Customization are executed in a positive sequence.
                6: [('salv', ([5, 6], 151.1256221)), ...],  # Set A Local Value whose local indices are (5,6) in element 6 to
                    151.1256221.
                ...
        }
        """
        return self.___customizations___

    def ___PRIVATE_do_execute_customization___(self, RETURN, e):
        """Execute all added customization.

        :param RETURN: The RETURN to be customized.
        :param e: for #e element.
        :return:
        """

        if e not in self._customizations_:

            pass

        else:
            assert spspa.isspmatrix_csc(RETURN) or spspa.isspmatrix_csr(RETURN)

            CUSi = self._customizations_[e] # customizations for the output of #e element.


            for cus in CUSi:

                key, factors = cus # for example cus=('clr', 5) or ('salv', ([1,1],5))

                # clear a local row of the EWC-sparse-matrix ----------------- BELOW -----------------------------------
                if key == 'clr': # Clear Local Row
                    # `factors` will be an int
                    assert factors % 1 == 0, f"factors={factors} wrong, when `clr`, factors must be an int."
                    if not spspa.isspmatrix_lil(RETURN): RETURN = RETURN.tolil()
                    RETURN[factors, :] = 0

                # clear local rows for element #e --------------------- BELOW ---------------------------
                elif key == 'clrs': # Clear Local Row
                    #factores are rows
                    if not spspa.isspmatrix_lil(RETURN): RETURN = RETURN.tolil()
                    RETURN[factors, :] = 0

                # Set A Local Value of local indices (i, j) to v --------------------- BELOW ---------------------------
                elif key == 'salv':
                    ij, v = factors
                    i, j = ij
                    if not spspa.isspmatrix_lil(RETURN): RETURN = RETURN.tolil()
                    RETURN[i, j] = v

                # identify Local Rows --------------------- BELOW ---------------------------
                elif key == 'ilrs':
                    if not spspa.isspmatrix_lil( RETURN): RETURN = RETURN.tolil()
                    RETURN[factors, :] = 0
                    RETURN[factors, factors] = 1

                # identify Local Rows At Columns --------------------- BELOW ---------------------------
                elif key == 'ilrsac':
                    # key, factors = ('ilrsac', ([0,1,2,3], [10, 11, 12 ,14]))
                    if not spspa.isspmatrix_lil( RETURN): RETURN = RETURN.tolil()
                    RETURN[factors[0], :] = 0
                    RETURN[factors[0], factors[1]] = 1

                # Not Implemented --------------------------------------------- BELOW ----------------------------------
                else:
                    raise NotImplementedError(f"Can not handle customization key={key}.")
                #============================================================== ABOVE ==================================


        if not (spspa.isspmatrix_csc(RETURN) or spspa.isspmatrix_csr(RETURN)): RETURN = RETURN.tocsr()


        return RETURN




    def clear_global_row(self, r):
        """Make the #r row of the global matrix (assemble self to get the global matrix) to be all zero.
        To achieve this, we only adjust the local output of the call function.

        We need to first find which elements contribute to global row #r, and then find the "local rows".
        Afterwards, we can apply customization key 'clr' with customization factor "local row" to elements.

        :param r:
        :return:
        """
        assert not self._spa_mat_.do.___locker___, f"the assembled matrix is locked!"
        assert not self._spa_mat_.do.___sparsity_locker___, f"the sparsity is locked!"

        assert r % 1 == 0, f"r={r}({r.__class__.__name__}) is wrong."
        if isinstance(r, float): r = int(r)

        rCGM = self._spa_mat_.gathering_matrices[0]
        assert rCGM is not None, "I have no row gathering_matrix!"
        OUTPUT = rCGM.do.find.elements_and_local_indices_of_dof(r)
        if OUTPUT is None:
            pass
        else:
            elements, local_indices = OUTPUT
            for e, i in zip(elements, local_indices):

                assert e in self._spa_mat_, f"element {e} is not in this core!"

                if e not in self._customizations_:
                    self._customizations_[e] = list()

                self._customizations_[e].append(('clr', int(i)))

    def clear_global_rows(self, rs):
        """

        :param rs:
        :return:
        """
        assert not self._spa_mat_.do.___locker___, f"the assembled matrix is locked!"
        assert not self._spa_mat_.do.___sparsity_locker___, f"the sparsity is locked!"
        raise NotImplementedError(f"{rs}")



    def set_assembled_M_ij_to(self, i, j, v):
        """Let M be the assembled matrix, we set M[i,j] = v.

        :param i:
        :param j:
        :param v:
        :return:
        """
        assert not self._spa_mat_.do.___locker___, f"the assembled matrix is locked!"
        assert not self._spa_mat_.do.___sparsity_locker___, f"the sparsity is locked!"

        assert i % 1 == 0, f"i={i}({i.__class__.__name__}) is wrong."
        if isinstance(i, float): i = int(i)
        assert j % 1 == 0, f"j={j}({j.__class__.__name__}) is wrong."
        if isinstance(j, float): j = int(j)

        rCGM, cCGM = self._spa_mat_.gathering_matrices
        assert rCGM is not None, "I have no row gathering_matrix!"
        assert cCGM is not None, "I have no col gathering_matrix!"

        rO = rCGM.do.find.elements_and_local_indices_of_dof(i)
        cO = cCGM.do.find.elements_and_local_indices_of_dof(j)

        if rO is not None and cO is not None:
            rE, rI = rO
            cE, cI = cO

            TBG = list()

            LOCAL_ij = list()

            for i, re in enumerate(rE):
                for j, ce in enumerate(cE):
                    if ce == re:
                        TBG.append(re)
                        LOCAL_ij.append((rI[i], cI[j]))
                        break
                    else:
                        pass

            if TBG == list():
                TBG = None

        else:
            TBG = None # to be gathered to master

        ALL_TBG = COMM.gather(TBG, root=MASTER_RANK)
        if RANK == MASTER_RANK:
            ELEMENTS_HAVE_THE_TARGET = list()
            for tbg in ALL_TBG:
                if tbg is not None:
                    ELEMENTS_HAVE_THE_TARGET.extend(tbg)
                else:
                    pass
            ELEMENTS_HAVE_THE_TARGET = list(set(ELEMENTS_HAVE_THE_TARGET))

            assert len(ELEMENTS_HAVE_THE_TARGET) > 0, "We are trying to set value to an empty place at which the " \
                                                      "assembling will allocate no value."


            the_element = ELEMENTS_HAVE_THE_TARGET[0] # we must find an element by here.
        else:
            the_element = None
        the_element = COMM.bcast(the_element, root=MASTER_RANK)


        if TBG is not None:
            if the_element in TBG:
                I_am_in = RANK
            else:
                I_am_in = SIZE
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
            assert the_element in TBG, "Must be like this!"
            for k, e in enumerate(TBG):
                # noinspection PyUnboundLocalVariable
                i, j = LOCAL_ij[k]
                if e not in self.___customizations___:
                    self.___customizations___[e] = list()
                if e != the_element:
                    self.___customizations___[e].append(('salv', ([i, j], 0)))
                else:
                    self.___customizations___[e].append(('salv', ([i, j], v)))


        else: # These cores have no business with the setting the value, they only need to make the value to be zero.
            if TBG is None: # these cores have no business at all.
                pass
            else:
                for k, e in enumerate(TBG):
                    # noinspection PyUnboundLocalVariable
                    i, j = LOCAL_ij[k]
                    if e not in self.___customizations___:
                        self.___customizations___[e] = list()
                    self.___customizations___[e].append(('salv', ([i, j], 0)))




    def identify_global_row(self, r):
        """Let M be the assembled matrix. We want to make M[r,:] = 0 except M[r,r] = 1.

        :param r:
        :return:
        """
        assert not self._spa_mat_.do.___locker___, f"the assembled matrix is locked!"
        assert not self._spa_mat_.do.___sparsity_locker___, f"the sparsity is locked!"

        assert r % 1 == 0, f"r={r}({r.__class__.__name__}) is wrong."
        if isinstance(r, float): r = int(r)
        self.clear_global_row(r)
        self.set_assembled_M_ij_to(r, r, 1)

    def identify_global_rows(self, rs):
        """

        :param rs:
        :return:
        """
        assert not self._spa_mat_.do.___locker___, f"the assembled matrix is locked!"
        assert not self._spa_mat_.do.___sparsity_locker___, f"the sparsity is locked!"

        raise NotImplementedError(f'{rs}')



    def identify_global_rows_according_to(self, i, interpret, AS='local'):
        """

        Parameters
        ----------
        i
        interpret :
            for example, `u.BC.interpret`.
        AS

        Returns
        -------

        """
        assert not self._spa_mat_.do.___locker___, f"the assembled matrix is locked!"
        assert not self._spa_mat_.do.___sparsity_locker___, f"the sparsity is locked!"

        assert self._spa_mat_.elements._mesh_ == interpret._mesh_, "BC_itp mesh != EWC mesh."
        bsp = self._spa_mat_.bmat_shape

        assert bsp is not False and np.shape(bsp) == (2,), \
            "we can only use this to bmat EWC matrices to keep the block safe. " \
            "if you want to apply it to a single EWC matrix, first wrap it in a list, then send it to bmat"

        assert i % 1 == 0, f"i={i}({i.__class__.__name__}) cannot be an index!"
        if isinstance(i, float): i = int(i)
        assert 0 <= i < bsp[0], f"i={i} is beyond the shape of this " \
                                f"EWC spare matrix of block shape ({bsp[0]},.)."

        assert self._spa_mat_._gathering_matrices_0_ is not None, \
            f"We at least need row gathering matrix!"

        if AS == 'local':
            LDF = interpret.local.dofs

            GM = self._spa_mat_.gathering_matrices[0]
            if GM.___Pr_IS_regular___:
                LDR = GM.___Pr_regular__local_dofs_ranges___[i]
                start = LDR.start
            else:
                raise NotImplementedError(f"Not implemented for irregular Gathering Matrix.")

            for e in LDF: # go through all locally involved mesh element numbers
                assert e in self._spa_mat_

                if e not in self.___customizations___:
                    self.___customizations___[e] = list()

                local_dofs = np.array(LDF[e], dtype=int) + start

                self.___customizations___[e].append(('ilrs', local_dofs))

        else:
            raise Exception(f"Cannot identify global rows through"
                            f" interpret interpreted"
                            f" as <{AS}>.")

    def off_diagonally_identify_rows_according_to(
        self, i, j, row_itp, col_itp, AS='local'):
        """We will identify off-diagonal block, block[i][j], and set
        block[i][:][row_pds.dofs, :] = 0
        block[i][j][row_pds.dofs, col_pds.dofs] = 1
        locally. Note it is "locally", that means, after assembly, we
        may have values greater than 1. This is okay because same thing
        will happen in the right-hand-side EWC vector.

        :param i:
        :param j:
        :param row_itp:
        :param col_itp:
        :param AS:
        :return:
        """
        assert not self._spa_mat_.do.___locker___, f"the assembled matrix is locked!"
        assert not self._spa_mat_.do.___sparsity_locker___, f"the sparsity is locked!"

        assert self._spa_mat_.elements._mesh_ == row_itp._mesh_, \
            "row PartialDofs elements != EWC elements."
        assert self._spa_mat_.elements._mesh_ == col_itp._mesh_, \
            "col PartialDofs elements != EWC elements."


        assert self._spa_mat_._gathering_matrices_0_ is not None, \
            f"We need row gathering matrix!"
        assert self._spa_mat_._gathering_matrices_1_ is not None, \
            f"We need col gathering matrix!"

        bsp = self._spa_mat_.bmat_shape
        assert bsp is not False and np.shape(bsp) == (2,), \
            "we can only use this to bmat EWC matrices to keep the block safe. " \
            "if you want to apply it to a single EWC matrix, first wrap it in a list, then send it to bmat"


        assert i % 1 == 0, f"i={i}({i.__class__.__name__}) cannot be an index!"
        if isinstance(i, float): i = int(i)
        assert 0 <= i < bsp[0], f"i={i} is beyond the shape of this " \
                                f"EWC spare matrix of block shape ({bsp[0]},.)."

        assert j % 1 == 0, f"j={j}({j.__class__.__name__}) cannot be an index!"
        if isinstance(j, float): j = int(j)
        assert 0 <= j < bsp[1], f"j={j} is beyond the shape of this " \
                                f"EWC spare matrix of block shape ({bsp[1]},.)."

        if AS == 'local':

            GM = self._spa_mat_.gathering_matrices
            if GM[0].___Pr_IS_regular___:
                LDR_row = GM[0].___Pr_regular__local_dofs_ranges___[i]
                start_row = LDR_row.start
            else:
                raise NotImplementedError(f"Not implemented for irregular Gathering Matrix.")

            if GM[1].___Pr_IS_regular___:
                LDR_col = GM[1].___Pr_regular__local_dofs_ranges___[j]
                start_col = LDR_col.start
            else:
                raise NotImplementedError(f"Not implemented for irregular Gathering Matrix.")

            LDF_row = row_itp.local.dofs
            LDF_col = col_itp.local.dofs

            for e in LDF_row: # go through all locally involved mesh element numbers
                assert e in self._spa_mat_ and e in LDF_col

                if e not in self.___customizations___:
                    self.___customizations___[e] = list()

                ROW = LDF_row[e]
                COL = LDF_col[e]

                row_local_dofs = np.array(ROW, dtype=int) + start_row
                col_local_dofs = np.array(COL, dtype=int) + start_col

                self.___customizations___[e].append(
                    ('ilrsac', (row_local_dofs, col_local_dofs)))

        else:
            raise Exception(f"Cannot off-diagonally identify global rows through "
                            f"interpret interpreted "
                            f"as <{AS}>.")
