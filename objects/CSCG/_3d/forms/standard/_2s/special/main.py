# -*- coding: utf-8 -*-


import sys
if './' not in sys.path: sys.path.append('./')

from objects.CSCG._3d.forms.standard._2s.special.vortex_detection import ___3dCSCG_2Form_Vortex_Detection___

from screws.freeze.main import FrozenOnly
from tools.linear_algebra.elementwise_cache.objects.sparse_matrix.main import EWC_SparseMatrix
from tools.linear_algebra.elementwise_cache.objects.column_vector.main import EWC_ColumnVector

from objects.CSCG._3d.forms.standard._2s.special.helpers.cross_product_2__ip_2_2M0 import ___3dCSCG_2Form_CrossProduct_2__ip_2_2M0___
from objects.CSCG._3d.forms.standard._2s.special.helpers.cross_product_2__ip_2 import ___3dCSCG_2Form_CrossProduct_2__ip_2___
from objects.CSCG._3d.forms.standard._2s.special.helpers.cross_product_2__ip_1 import ___3dCSCG_2Form_CrossProduct_2__ip_1___




class _2Form_Special(FrozenOnly):
    def __init__(self, _2sf):
        self._sf_ = _2sf
        self._vortex_detection_ = None
        self._freeze_self_()

    def cross_product_2f__ip_2f(self, u, e, quad_degree=None, output='2-M-1', cache=None):
        """
        (self X 2form, 2form)

        (self X u    , e    )

        We do ``(self X other, e)`` where ``self`` and ``other`` both are n-form, n be either 1 or 2.

        output:
            '2-M-1': Means we return a local matrix refers to local dofs of e (rows, index-0) and u (cols, index-1)
            'MDM' : return a Multi-Dimensional-Matrix whose 0-, 1-, 2-axis refer to self, u, e.

        cache :
            A cache to pass known variables (to save computational time and memory cost.)

        :return:
        """
        if output == '2-M-1':
            SCP_generator = ___3dCSCG_2Form_CrossProduct_2__ip_2___(self._sf_, u, e, quad_degree=quad_degree)
        elif output == 'MDM':
            SCP_generator = ___3dCSCG_2Form_CrossProduct_2__ip_2___(self._sf_, u, e, quad_degree=quad_degree)
            return SCP_generator.MDM
        elif output == '2-M-0':
            SCP_generator = ___3dCSCG_2Form_CrossProduct_2__ip_2_2M0___(self._sf_, u, e, quad_degree=quad_degree,
                                                                        cache=cache)
        else:
            raise NotImplementedError(f"output={output} is not implemented.")

        return EWC_SparseMatrix(self._sf_.mesh.elements, SCP_generator, 'no_cache')

    def cross_product_2f__ip_1f(self, u, e, quad_degree=None, output='2-M-1'):
        """
        (self X 2form, 1form)

        (self X u    , e    )

        We do ``(self X other, e)`` where ``self`` and ``other`` both are 2-form.

        output:
            '2-M-1': Means we return a local matrix refers to local dofs of e (rows, index-0) and u (cols, index-1)

        :return:
        """
        if output == '2-M-1':
            SCP_generator = ___3dCSCG_2Form_CrossProduct_2__ip_1___(self._sf_, u, e, quad_degree=quad_degree)
        else:
            raise NotImplementedError(f"output={output} is not implemented.")

        return EWC_SparseMatrix(self._sf_.mesh.elements, SCP_generator, 'no_cache')

    @property
    def vortex_detection(self):
        if self._vortex_detection_ is None:
            self._vortex_detection_ = ___3dCSCG_2Form_Vortex_Detection___(self._sf_)
        return self._vortex_detection_

    def hybrid_pairing(self, adt2, time=0):
        """"""
        assert adt2.__class__.__name__ == '_3dCSCG_T2_ADF', f"I need a 3dCSCG AD-Trace-2-form."

        sf = self._sf_
        mesh = sf.mesh

        assert sf.TW.BC.body is not None, f'3dCSCG primal 2-sf has no TW.BC function.'
        assert sf.BC.valid_boundaries is not None, f'3dCSCG primal 2-sf has no valid boundary.'
        assert adt2.prime.TW.BC.body is not None, f'3dCSCG ad-2-trace has no TW.BC function.'
        assert adt2.BC.valid_boundaries is not None, f'3dCSCG ad-2-trace has no valid boundary.'

        sf.TW.do.push_BC_to_instant(time)
        adt2.prime.TW.do.push_BC_to_instant(time)

        T = adt2.matrices.trace
        D = EWC_SparseMatrix(mesh, (adt2.num.basis, adt2.num.basis))
        b = EWC_ColumnVector(mesh, adt2)

        T.gathering_matrices = (adt2, sf)
        D.gathering_matrices = (adt2, adt2)
        b.gathering_matrix = adt2

        #----- get boundaries and do a check --------------------------------------
        Dirichlet_boundaries = adt2.BC.valid_boundaries
        Neumann_boundaries = sf.BC.valid_boundaries

        bns = mesh.boundaries.names
        SDb = set(Dirichlet_boundaries)
        SNb = set(Neumann_boundaries)
        assert SDb & SNb == set()   , f"Dirichlet_boundaries intersect Neumann_boundaries is not None."
        assert SDb | SNb == set(bns), f"Dirichlet_boundaries union Neumann_boundaries is not full!"

        #-------- set Neumann boundary condition ---------------------------------------------------
        sf.BC.valid_boundaries = Neumann_boundaries
        adt2.BC.valid_boundaries = Neumann_boundaries
        col_pc = sf.BC.partial_cochain
        row_pd = adt2.BC.partial_dofs
        T = T.adjust.identify_rows_according_to_two_CSCG_partial_dofs(row_pd, col_pc)
        b = b.adjust.set_entries_according_to_CSCG_partial_cochains(row_pd, col_pc)

        #-------- set Dirichlet boundary condition -------------------------------
        adt2.BC.valid_boundaries = Dirichlet_boundaries
        adt_pc = adt2.BC.partial_cochain
        D = D.adjust.identify_rows_according_to_CSCG_partial_dofs(adt_pc)
        T = T.adjust.clear_rows_according_to_CSCG_partial_dofs(adt_pc)
        b = b.adjust.set_entries_according_to_CSCG_partial_cochains(adt_pc, adt_pc)

        #=====================================================================================
        return T, D, b