

from screws.freeze.main import FrozenOnly
from screws.decorators.accepts import accepts
from scipy import sparse as spspa
from root.config.main import *


class SpaMat_Customize(FrozenOnly):
    """Store all method to customize the EWC sparse matrix."""
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


    @accepts('self', (int, float, 'int32', 'int64'))
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
        raise NotImplementedError()


    @accepts('self',
             (int, float, 'int32', 'int64'),
             (int, float, 'int32', 'int64'),
             (int, float, 'int32', 'int64'))
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

        ALL_TBG = cOmm.gather(TBG, root=mAster_rank)
        if rAnk == mAster_rank:
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
        the_element = cOmm.bcast(the_element, root=mAster_rank)


        if TBG is not None:
            if the_element in TBG:
                I_am_in = rAnk
            else:
                I_am_in = sIze
        else:
            I_am_in = sIze

        who_are_in = cOmm.gather(I_am_in, root=mAster_rank)
        if rAnk == mAster_rank:
            use_who = min(who_are_in)
            assert use_who != sIze, "Something is wrong, we need to find a correct core."
        else:
            use_who = None
        use_who = cOmm.bcast(use_who, root=mAster_rank)

        if use_who == rAnk:
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


    @accepts('self', (int, float, 'int32', 'int64'))
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

        raise NotImplementedError()


    @accepts('self',
             (int, float, 'int32', 'int64'),
             ('PartialDofs', 'PartialCochain'))
    def identify_global_rows_according_to_CSCG_partial_dofs(self, i, pds, interpreted_as='local_dofs'):
        """We locally make M[r,:] = 0 except M[r,r] = 1 according to a CSCG PartialDofs
        instance in block[i][:]

        This will ask the local and global matrix to be square: But we
        do not check it.

        IMPORTANT: as this is done only locally, after assembly, M[r,r] may be greater than 1.
            for example, for a corner point, M[r,r] can be 8 if it is internal. This, for example,
            will not introduce problem for imposing bc if in the right-hand-side vector, we
            also give a 8 times greater value.

        :param pds: we will get dofs form this.
        :param i: We locally make block[i][:][dofs, :] = 0 except block[i][i][dofs, dofs] = 1
        :param interpreted_as:
        :return:
        """
        assert not self._spa_mat_.do.___locker___, f"the assembled matrix is locked!"
        assert not self._spa_mat_.do.___sparsity_locker___, f"the sparsity is locked!"

        if pds.__class__.__name__ == 'PartialCochain': pds = pds.dofs

        assert self._spa_mat_.elements._mesh_ == pds._mesh_, \
            "PartialDofs elements != EWC elements."

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


        if interpreted_as == 'local_dofs':

            GM = self._spa_mat_.gathering_matrices[0]
            LDR = GM.local_dofs_ranges[i]
            start = LDR.start

            LDF = pds.interpreted_as.local_dofs

            for e in LDF: # go through all locally involved mesh element numbers
                assert e in self._spa_mat_

                if e not in self.___customizations___:
                    self.___customizations___[e] = list()

                local_dofs = np.array(LDF[e]) + start

                self.___customizations___[e].append(('ilrs', local_dofs))

        else:
            raise Exception(f"Cannot identify global rows through "
                            f"SCG_partial_dofs interpreted "
                            f"as <{interpreted_as}>.")



    @accepts('self',
             (int, float, 'int32', 'int64'),
             (int, float, 'int32', 'int64'),
             ('PartialDofs', 'PartialCochain'),
             ('PartialDofs', 'PartialCochain'))
    def off_diagonally_identify_rows_according_to_two_CSCG_partial_dofs(
        self, i, j, row_pds, col_pds, interpreted_as='local_dofs'):
        """We will identify off-diagonal block, block[i][j], and set
        block[i][:][row_pds.dofs, :] = 0
        block[i][j][row_pds.dofs, col_pds.dofs] = 1
        locally. Note it is "locally", that means, after assembly, we
        may have values greater than 1. This is okay because same thing
        will happen in the right-hand-side EWC vector.

        :param i:
        :param j:
        :param row_pds:
        :param col_pds:
        :param interpreted_as:
        :return:
        """
        assert not self._spa_mat_.do.___locker___, f"the assembled matrix is locked!"
        assert not self._spa_mat_.do.___sparsity_locker___, f"the sparsity is locked!"

        if row_pds.__class__.__name__ == 'PartialCochain':
            row_pds = row_pds.dofs
        if col_pds.__class__.__name__ == 'PartialCochain':
            col_pds = col_pds.dofs

        assert self._spa_mat_.elements._mesh_ == row_pds._mesh_, \
            "row PartialDofs elements != EWC elements."
        assert self._spa_mat_.elements._mesh_ == col_pds._mesh_, \
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

        if interpreted_as == 'local_dofs':

            GM = self._spa_mat_.gathering_matrices
            LDR_row = GM[0].local_dofs_ranges[i]
            LDR_col = GM[1].local_dofs_ranges[j]
            start_row = LDR_row.start
            start_col = LDR_col.start

            LDF_row = row_pds.interpreted_as.local_dofs
            LDF_col = col_pds.interpreted_as.local_dofs

            # print(GM[0].local_dofs_ranges, start_row, start_col)

            for e in LDF_row: # go through all locally involved mesh element numbers
                assert e in self._spa_mat_ and e in LDF_col

                if e not in self.___customizations___:
                    self.___customizations___[e] = list()

                ROW = LDF_row[e]
                COL = LDF_col[e]

                COL = self.___PRIVATE_correcting_correspondence___(
                    row_pds._form_, ROW, COL, col_pds._form_)

                row_local_dofs = np.array(ROW) + start_row
                col_local_dofs = np.array(COL) + start_col

                self.___customizations___[e].append(
                    ('ilrsac', (row_local_dofs, col_local_dofs)))

        else:
            raise Exception(f"Cannot off-diagonally identify global rows through "
                            f"SCG_partial_dofs interpreted "
                            f"as <{interpreted_as}>.")

    def ___PRIVATE_correcting_correspondence___(self, rf, R, C, cf):
        """

        Parameters
        ----------
        rf
        R
        C
        cf

        Returns
        -------

        """
        if rf.__class__.__name__ == '_3dCSCG_1Trace' and cf.__class__.__name__ == '_3dCSCG_1Form':
            # we have to make sure this to make the singularity handling possible
            if self._t1f_s1f_cc_ is None:
                s1f = list()
                for side in 'NSWEBF':
                    s1f.append(cf.numbering.do.find.local_dofs_on_element_side(side))
                s1f = np.concatenate(s1f)

                self._t1f_s1f_cc_ = s1f
            else:
                s1f = self._t1f_s1f_cc_

            cC = s1f[R]
            assert len(cC) == len(C) and set(cC) == set(C), f"must be this case."
            return cC

        elif rf.__class__.__name__ == '_3dCSCG_0Trace' and cf.__class__.__name__ == '_3dCSCG_0Form':
            # we have to make sure this to make the singularity handling possible
            if self._t0f_s0f_cc_ is None:
                s0f = list()
                for side in 'NSWEBF':
                    s0f.append(cf.numbering.do.find.local_dofs_on_element_side(side))
                s0f = np.concatenate(s0f)

                self._t0f_s0f_cc_ = s0f
            else:
                s0f = self._t0f_s0f_cc_

            cC = s0f[R]
            assert len(cC) == len(C) and set(cC) == set(C), f"must be this case."
            return cC

        else:
            return C
