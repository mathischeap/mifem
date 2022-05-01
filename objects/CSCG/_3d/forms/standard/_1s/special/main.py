
import sys
if './' not in sys.path: sys.path.append('./')

from numpy import array
from itertools import chain
from screws.freeze.main import FrozenOnly
from tools.linear_algebra.gathering.regular.matrix.main import Gathering_Matrix
from tools.linear_algebra.gathering.vector import Gathering_Vector

from tools.linear_algebra.elementwise_cache.objects.sparse_matrix.main import EWC_SparseMatrix
from tools.linear_algebra.elementwise_cache.objects.column_vector.main import EWC_ColumnVector

from objects.CSCG._3d.forms.standard._1s.special.vortex_detection import \
    ___3dCSCG_1Form_Vortex_Detection___
from objects.CSCG._3d.forms.standard._1s.special.helpers.cross_product_1__ip_1 import \
    ___3dCSCG_1Form_CrossProduct_1__ip_1___
from objects.CSCG._3d.forms.standard._1s.special.helpers.cross_product_2__ip_2 import \
    ___3dCSCG_1Form_CrossProduct_2__ip_2___
from root.config.main import cOmm, MPI




class _1Form_Special(FrozenOnly):
    def __init__(self, _1sf):
        self._sf_ = _1sf
        self._vortex_detection_ = None
        self._freeze_self_()

    def cross_product_1f__ip_1f(self, u, e, quad_degree=None, output='2-M-1'):
        """
        (self X 1-form, 1-form). To first cross product with a 1-form then do an inner product with
        another 1-form.

        output:
            '2-M-1': Means we return a local matrix refers to local dofs of e(column) and u (row)

        :return:
        """
        if output == '2-M-1':
            SCP_generator = ___3dCSCG_1Form_CrossProduct_1__ip_1___(self._sf_, u, e, quad_degree=quad_degree)
        else:
            raise NotImplementedError(f"output={output} is not implemented.")

        return EWC_SparseMatrix(self._sf_.mesh.elements, SCP_generator, 'no_cache')

    def cross_product_2f__ip_2f(self, u, e, quad_degree=None, output='2-M-1'):
        """
        (self X 2-form, 2-form). To first cross product with a 2-form then do an inner product with
        another 2-form.

        output:
            '2-M-1': Means we return a local matrix refers to local dofs of e (column) and u (row)

        :return:
        """
        if output == '2-M-1':
            SCP_generator = ___3dCSCG_1Form_CrossProduct_2__ip_2___(self._sf_, u, e, quad_degree=quad_degree)
        else:
            raise NotImplementedError(f"output={output} is not implemented.")

        return EWC_SparseMatrix(self._sf_.mesh.elements, SCP_generator, 'no_cache')

    @property
    def vortex_detection(self):
        if self._vortex_detection_ is None:
            self._vortex_detection_ = ___3dCSCG_1Form_Vortex_Detection___(self._sf_)
        return self._vortex_detection_

    def hybrid_pairing(self, adt1, e1, time=0):
        """"""
        assert adt1.__class__.__name__ == '_3dCSCG_T1_ADF', f"I need a 3dCSCG AD-Trace-1-form."
        assert e1.__class__.__name__ == '_3dCSCG_1Edge', f"I need a 3dCSCG 1-edge-form."
        sf = self._sf_
        mesh = sf.mesh

        assert sf.TW.BC.body is not None, f'3dCSCG primal 1-sf has no TW.BC function.'
        assert sf.BC.valid_boundaries is not None, f'3dCSCG primal 1-sf has no valid boundary.'
        assert adt1.prime.TW.BC.body is not None, f'3dCSCG ad-1-trace has no TW.BC function.'
        assert adt1.BC.valid_boundaries is not None, f'3dCSCG ad-1-trace has no valid boundary.'

        sf.TW.do.push_BC_to_instant(time)
        adt1.prime.TW.do.push_BC_to_instant(time)

        T = adt1.matrices.trace
        D = EWC_SparseMatrix(mesh, (adt1.num.basis, adt1.num.basis))
        C = e1.matrices.complement

        b = EWC_ColumnVector(mesh, adt1)

        T.gathering_matrices = (adt1, sf)
        D.gathering_matrices = (adt1, adt1)
        C.gathering_matrices = (adt1, e1)
        b.gathering_matrix = adt1

        #----- get boundaries and do a check --------------------------------------
        Dirichlet_boundaries = adt1.BC.valid_boundaries
        Neumann_boundaries = sf.BC.valid_boundaries

        bns = mesh.boundaries.names
        SDb = set(Dirichlet_boundaries)
        SNb = set(Neumann_boundaries)
        assert SDb & SNb == set()   , f"Dirichlet_boundaries intersect Neumann_boundaries is not None."
        assert SDb | SNb == set(bns), f"Dirichlet_boundaries union Neumann_boundaries is not full!"
        #-------- set Neumann boundary condition ---------------------------------------------------

        sf.BC.valid_boundaries = Neumann_boundaries
        adt1.BC.valid_boundaries = Neumann_boundaries
        col_pc = sf.BC.partial_cochain
        row_pd = adt1.BC.partial_dofs
        T = T.adjust.identify_rows_according_to_two_CSCG_partial_dofs(row_pd, col_pc)
        b = b.adjust.set_entries_according_to_CSCG_partial_cochains(row_pd, col_pc)

        #-------- set Dirichlet boundary condition -------------------------------
        adt1.BC.valid_boundaries = Dirichlet_boundaries
        adt_pc = adt1.BC.partial_cochain
        D = D.adjust.identify_rows_according_to_CSCG_partial_dofs(adt_pc)
        T = T.adjust.clear_rows_according_to_CSCG_partial_dofs(adt_pc)
        b = b.adjust.set_entries_according_to_CSCG_partial_cochains(adt_pc, adt_pc)

        #---------------- Send T, C for hybrid singularity overcoming ------------------------------
        T, C, SKIPPED_edge_elements = self.___PRIVATE_overcoming_hybrid_singularity___(
            T, C, Dirichlet_boundaries=Dirichlet_boundaries)

        #------------- make a special Gathering matrix for the 1-edge-form ------------------------
        eGM = self.___PRIVATE_1ef_hybrid_GM___(SKIPPED_edge_elements)

        return T, D, C, b, eGM

    def ___PRIVATE_1ef_hybrid_GM___(self, SKIPPED_edge_elements):
        """We make a special gathering matrix of the 1-edge-form using for the hybrid singularity overcoming."""
        mesh = self._sf_.mesh
        p = self._sf_.space.p
        px, py, pz = p

        D_p_D = {'NS':px, 'WE':py, 'BF':pz} # Direction-p-Dict

        #------- take care SKIPPED_edge_elements -----------------------------------------------
        SKD = dict()
        for e in SKIPPED_edge_elements:
            if e in mesh.edge.elements:
                meee = mesh.edge.elements[e]
                direction = meee.direction
                SKD[e] = D_p_D[direction]

        ___ = cOmm.allgather(SKD)
        SKD = dict()
        for _ in ___: SKD.update(_)

        #---------------------------------------------------------------------------------------
        _EEW_GV_ = dict()

        TA_NB = mesh.edge.elements.___PRIVATE_find_type_and_amount_numbered_before___()
        assert len(TA_NB) == len(mesh.edge.elements)

        MINUS = 0

        for ee in mesh.edge.elements:

            if ee in SKD:

                _EEW_GV_[ee] = [0 for _ in range(SKD[ee])]

            else:
                ta_NB = TA_NB[ee]
                nx, ny, nz = ta_NB
                BEFORE = nx * px + ny * py + nz * pz # normal numbering condition

                _2bd_ = list()
                for sk in SKD:
                    if sk < ee:
                        MINUS += SKD[sk]
                        _2bd_.append(sk)
                    else:
                        break

                for _ in _2bd_:
                    del SKD[_]

                BEFORE -= MINUS

                direction = mesh.edge.elements[ee].direction
                if direction == 'NS':
                    _EEW_GV_[ee] = range(BEFORE, BEFORE + px)
                elif direction == 'WE':
                    _EEW_GV_[ee] = range(BEFORE, BEFORE + py)
                elif direction == 'BF':
                    _EEW_GV_[ee] = range(BEFORE, BEFORE + pz)
                else:
                    raise Exception()

        eGM = dict()
        for me in mesh.elements:
            E_MAP = mesh.edge.elements.map[me]
            GV = list()
            for mee in E_MAP:
                GV.append(_EEW_GV_[mee])
            GV = array([_ for _ in chain(*GV)])

            eGM[me] = Gathering_Vector(me, GV)

        eGM = Gathering_Matrix(eGM, mesh_type='_3dCSCG')

        return eGM

    def ___PRIVATE_overcoming_hybrid_singularity___(self, T, C, Dirichlet_boundaries=None):
        """

        Parameters
        ----------
        T :
            The trace matrix.
        C :
            The complementary matrix.
        Dirichlet_boundaries :
            The mesh boundaries where we will apply direct boundary condition to the (AD)trace dofs.
            For example, in the Poisson problem. The potential boundaries are a Dirichlet_boundaries.

        Returns
        -------

        """
        assert self._sf_.IS.hybrid, f"Only hybrid 1-form has this problem."
        assert T.__class__.__name__ == 'EWC_SparseMatrix'
        assert C.__class__.__name__ == 'EWC_SparseMatrix'

        mesh = self._sf_.mesh

        boundaries_names = mesh.boundaries.names
        if Dirichlet_boundaries is None:
            Dirichlet_boundaries = list()
        elif isinstance(Dirichlet_boundaries, str):
            Dirichlet_boundaries = [Dirichlet_boundaries,]
        else:
            assert isinstance(Dirichlet_boundaries, (list, tuple)), \
                f"Dirichlet_boundaries={Dirichlet_boundaries} is wrong."
            assert len(set(Dirichlet_boundaries)) == len(Dirichlet_boundaries), \
                f"Repeated boundaries found in {Dirichlet_boundaries}."

        for _, Db in enumerate(Dirichlet_boundaries):
            assert Db in boundaries_names, \
                f"Dirichlet_boundaries[{_}] = {Dirichlet_boundaries[_]} is not a valid boundary."

        nT = dict() # the new Trace matrix.
        nC = dict() # the new Trace matrix.

        SKIPPED_edge_elements = list()

        for i in range(mesh.edge.elements.GLOBAL_num):

            if i in mesh.edge.elements:
                edge_element = mesh.edge.elements[i]
                on_mesh_boundaries = edge_element.on_mesh_boundaries
                skip = False
                for mb in on_mesh_boundaries:
                    if mb in Dirichlet_boundaries:
                        skip = True
                        break
            else:
                skip = False

            skip = cOmm.allreduce(skip, op=MPI.LOR)

            if skip:

                SKIPPED_edge_elements.append(i)

            else:

                SOS = mesh.edge.elements.do.find.hybrid_singularity_overcoming_setting(i)

                if SOS is None: # this core has no business with this SOS
                    pass
                else:
                    replacing = SOS.replacing
                    mesh_element, corner_edge = replacing

                    through = SOS.through
                    assert mesh_element in mesh.elements
                    sf_local_dofs = self._sf_.numbering.do.find.local_dofs_on_element_corner_edge(
                        corner_edge)

                    trace_element, trace_edge = through
                    T_MAP = mesh.trace.elements.map[mesh_element]
                    for si, _ in enumerate(T_MAP):
                        if _ == trace_element:
                            break
                    trace_face = 'NSWEBF'[si]
                    tf_local_dofs = self._sf_.space.local_numbering.\
                        ___PRIVATE_find_MESH_ELEMENT_WISE_local_dofs_of_1Trace_edge___(
                        trace_face, trace_edge
                    )

                    positions = edge_element.positions
                    for pos in positions:
                        if int(pos[:-2]) == mesh_element:
                            edge_name = pos[-2:]
                            break
                    ef_local_dofs = self._sf_.space.local_numbering.\
                        ___PRIVATE_find_MESH_ELEMENT_WISE_local_dofs_of_1edge_edge___(
                        edge_name
                    )

                    assert len(sf_local_dofs) == len(tf_local_dofs) == len(ef_local_dofs), \
                        f"Trivial check!"

                    if mesh_element not in nT:
                        nT[mesh_element] = T[mesh_element].copy().tolil()

                    V = nT[mesh_element][tf_local_dofs, sf_local_dofs]
                    nT[mesh_element][tf_local_dofs, sf_local_dofs] = 0

                    if mesh_element not in nC:
                        nC[mesh_element] = C[mesh_element].copy().tolil()
                    nC[mesh_element][tf_local_dofs, ef_local_dofs] = V

        for _ in T:
            if _ not in nT:
                nT[_] = T[_]
            else:
                # noinspection PyUnresolvedReferences
                nT[_] = nT[_].tocsr()

        for _ in C:
            if _ not in nC:
                nC[_] = C[_]
            else:
                # noinspection PyUnresolvedReferences
                nC[_] = nC[_].tocsr()

        nT = T.__class__(mesh, nT, cache_key_generator = 'no_cache')
        nC = C.__class__(mesh, nC, cache_key_generator = 'no_cache')

        return nT, nC, SKIPPED_edge_elements









if __name__ == '__main__':
    # mpiexec -n 5 python objects\CSCG\_3d\forms\standard\_1s\special\main.py
    from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller, ExactSolutionSelector

    elements = [2,2,2]
    # mesh = MeshGenerator('crazy_periodic', c=0.1)(elements)
    mesh = MeshGenerator('crazy', c=0.1)(elements)
    ES = ExactSolutionSelector(mesh)('Poisson:sincos1')

    Dirichlet_boundaries = ['Back', 'Front', 'West', ] #
    Neumann_boundaries = ["East", 'South', 'North', ]
    # mesh = MeshGenerator('bridge_arch_cracked')(elements)
    space = SpaceInvoker('polynomials')([2, 3, 4])
    FC = FormCaller(mesh, space)

    f1 = FC('1-f', is_hybrid=True)
    t1 = FC('1-adt')
    e1 = FC('1-e')

    f1.TW.BC.body = ES.status.velocity
    f1.TW.do.push_BC_to_instant(0)
    f1.BC.valid_boundaries = Neumann_boundaries

    t1.prime.TW.BC.body = ES.status.velocity.components.T_perp
    t1.prime.TW.do.push_BC_to_instant(0)
    t1.BC.valid_boundaries = Dirichlet_boundaries

    T, D, C, b, eGM = f1.special.hybrid_pairing(t1, e1)


    # T = t1.matrices.trace
    # C = e1.matrices.complement
    # T, C = f1.special.___PRIVATE_overcoming_hybrid_singularity___(
    #     T, C, Dirichlet_boundaries=Dirichlet_boundaries)[:2]

    # # # f1.dofs.visualize.matplot.connection_through_trace_dof(55, T, C, t1, e1, checking_mode=True)
    # #
    # #
    # # #
    # for i in range(t1.prime.numbering.gathering.GLOBAL_num_dofs):
    #     f1.dofs.visualize.matplot.connection_through_trace_dof(i, T, C, t1, e1, checking_mode=True)
    #
    # for i in range(e1.numbering.gathering.GLOBAL_num_dofs):
    #     f1.dofs.visualize.matplot.connection_through_around_edge_dof(
    #         i, T, C, t1, e1, checking_mode=True)
