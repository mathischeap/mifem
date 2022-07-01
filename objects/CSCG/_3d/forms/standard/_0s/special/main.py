# -*- coding: utf-8 -*-
import sys
if './' not in sys.path: sys.path.append('./')

from tools.linear_algebra.elementwise_cache.objects.sparse_matrix.main import EWC_SparseMatrix
from tools.linear_algebra.elementwise_cache.objects.column_vector.main import EWC_ColumnVector

from tools.linear_algebra.gathering.regular.matrix.main import Gathering_Matrix
from tools.linear_algebra.gathering.vector import Gathering_Vector

from itertools import chain
from screws.freeze.main import FrozenOnly
from root.config.main import np, cOmm, MPI, rAnk, mAster_rank
from objects.CSCG._3d.mesh.node.elements.do.find.helpers.SOS_internal import _3dCSCG_InternalNodeSOS
from objects.CSCG._3d.mesh.node.elements.do.find.helpers.SOS_boundary.corner import _3dCSCG_CornerNodeSOS


# from tqdm import tqdm
# from time import sleep
# from objects.CSCG._3d.mesh.node.elements.do.find.helpers.SOS_boundary.corner_edge import _3dCSCG_CornerEdgeNodeSOS
# from objects.CSCG._3d.mesh.node.elements.do.find.helpers.SOS_boundary.surface_middle import _3dCSCG_SurfaceMiddleNodeSOS




class _0Form_Special(FrozenOnly):
    def __init__(self, _0sf):
        self._sf_ = _0sf
        self._freeze_self_()

    def hybrid_pairing(self, adt0, e0, time=0):
        """"""
        assert adt0.__class__.__name__ == '_3dCSCG_T0_ADF', f"I need a 3dCSCG AD-Trace-0-form."
        assert e0.__class__.__name__ == '_3dCSCG_0Edge', f"I need a 3dCSCG 0-edge-form."
        sf = self._sf_
        mesh = sf.mesh

        assert sf.TW.BC.body is not None, f'3dCSCG primal 0-sf has no TW.BC function.'
        assert sf.BC.valid_boundaries is not None, f'3dCSCG primal 0-sf has no valid boundary.'
        assert adt0.prime.TW.BC.body is not None, f'3dCSCG ad-0-trace has no TW.BC function.'
        assert adt0.BC.valid_boundaries is not None, f'3dCSCG ad-0-trace has no valid boundary.'

        sf.TW.do.push_BC_to_instant(time)
        adt0.prime.TW.do.push_BC_to_instant(time)

        T = adt0.matrices.trace
        D = EWC_SparseMatrix(mesh, (adt0.num.basis, adt0.num.basis))
        C = e0.matrices.complement

        b = EWC_ColumnVector(mesh, adt0)

        T.gathering_matrices = (adt0, sf)
        D.gathering_matrices = (adt0, adt0)
        C.gathering_matrices = (adt0, e0)
        b.gathering_matrix = adt0

        #----- get boundaries and do a check --------------------------------------
        Dirichlet_boundaries = adt0.BC.valid_boundaries
        Neumann_boundaries = sf.BC.valid_boundaries

        bns = mesh.boundaries.names
        SDb = set(Dirichlet_boundaries)
        SNb = set(Neumann_boundaries)
        assert SDb & SNb == set()   , f"Dirichlet_boundaries intersect Neumann_boundaries is not None."
        assert SDb | SNb == set(bns), f"Dirichlet_boundaries union Neumann_boundaries is not full!"
        #-------- set Neumann boundary condition ---------------------------------------------------

        sf.BC.valid_boundaries = Neumann_boundaries
        adt0.BC.valid_boundaries = Neumann_boundaries
        col_pc = sf.BC.partial_cochain
        row_pd = adt0.BC.partial_dofs
        T = T.adjust.identify_rows_according_to_two_CSCG_partial_dofs(row_pd, col_pc)
        b = b.adjust.set_entries_according_to_CSCG_partial_cochains(row_pd, col_pc)

        #-------- set Dirichlet boundary condition -------------------------------
        adt0.BC.valid_boundaries = Dirichlet_boundaries
        adt_pc = adt0.BC.partial_cochain
        T = T.adjust.clear_rows_according_to_CSCG_partial_dofs(adt_pc)
        D = D.adjust.identify_rows_according_to_CSCG_partial_dofs(adt_pc)
        b = b.adjust.set_entries_according_to_CSCG_partial_cochains(adt_pc, adt_pc)

        T, C = self._sf_.special.___PRIVATE_overcoming_hybrid_singularity___(
            T, D, C, adt0, e0, Dirichlet_boundaries=Dirichlet_boundaries)

        eGM = self.___PRIVATE_hybrid_numbering_e0___(C)

        return T, D, C, b, eGM



    def ___PRIVATE_overcoming_hybrid_singularity___(
        self, T, D, C, adt0, e0, Dirichlet_boundaries=None):
        """Overcoming the hybrid singularity by adjusting trace matrix and complementary matrix.

        Parameters
        ----------
        T :
        D :
        C :
        Dirichlet_boundaries :

        Returns
        -------

        """
        assert self._sf_.IS.hybrid, f"Only hybrid 1-form has this problem."
        assert T.__class__.__name__ == 'EWC_SparseMatrix'
        assert D.__class__.__name__ == 'EWC_SparseMatrix'
        assert C.__class__.__name__ == 'EWC_SparseMatrix'

        assert adt0.__class__.__name__ == '_3dCSCG_T0_ADF', f"I need a 3dCSCG AD-Trace-0-form."
        assert e0.__class__.__name__ == '_3dCSCG_0Edge', f"I need a 3dCSCG 0-edge-form."

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

        T, C, SKIPPED_edge_elements = self.___PRIVATE_overcoming_hybrid_singularity_EDGE___(
            T, C, Dirichlet_boundaries)

        # T, C, SKIP_e0_dof = self.___PRIVATE_overcoming_hybrid_singularity_INTERNAL_NODE___(
        #     T, D, C, adt0, e0, Dirichlet_boundaries)

        return T, C


    def ___PRIVATE_overcoming_hybrid_singularity_EDGE___(self, T, C, Dirichlet_boundaries):
        """Overcoming the hybrid singularity on edge-elements.

        Parameters
        ----------
        T :
        C :
        Dirichlet_boundaries :

        Returns
        -------

        """
        mesh = self._sf_.mesh

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
                        corner_edge
                    )

                    trace_element, trace_edge = through
                    T_MAP = mesh.trace.elements.map[mesh_element]
                    for si, _ in enumerate(T_MAP):
                        if _ == trace_element:
                            break
                    trace_face = 'NSWEBF'[si]
                    tf_local_dofs = self._sf_.space.local_numbering.\
                        ___PRIVATE_find_MESH_ELEMENT_WISE_local_dofs_of_0Trace_edge___(
                        trace_face, trace_edge
                    )

                    positions = edge_element.positions
                    for pos in positions:
                        if int(pos[:-2]) == mesh_element:
                            edge_name = pos[-2:]
                            break
                    ef_local_dofs = self._sf_.space.local_numbering.\
                        ___PRIVATE_find_MESH_ELEMENT_WISE_local_dofs_of_0edge_edge___(
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


    def ___PRIVATE_hybrid_numbering_e0___(self, C):
        """"""
        sf = self._sf_
        mesh = sf.mesh
        assert len(C) == len(mesh.elements)

        E_MAP = mesh.edge.elements.map
        px, py, pz = sf.space.p

        numbered_edge_dof = dict()
        for i in C:
            Ci = C[i].tocsc()

            indptr = np.diff(Ci.indptr)

            MAP = E_MAP[i]

            start = 0
            for j, edge in enumerate(MAP):

                if j < 4:
                    p_ = px + 1
                elif j < 8:
                    p_ = py + 1
                else:
                    p_ = pz + 1

                if edge not in numbered_edge_dof:
                    numbered_edge_dof[edge] = np.zeros(p_)

                numbered_edge_dof[edge] += indptr[start : start + p_]

                start += p_

            assert len(indptr) == start

        numbered_edge_dof = cOmm.gather(numbered_edge_dof, root=mAster_rank)

        if rAnk == mAster_rank:
            NUMBER_EDGE_ELEMENTS = numbered_edge_dof[0]
            for nd in numbered_edge_dof[1:]:
                for i in nd:
                    if i in NUMBER_EDGE_ELEMENTS:
                        NUMBER_EDGE_ELEMENTS[i] += nd[i]
                    else:
                        NUMBER_EDGE_ELEMENTS[i] = nd[i]

            del numbered_edge_dof

            #--- start numbering --------------------------------------------------
            start = 0
            for i in range(mesh.edge.elements.GLOBAL_num):
                Ni = NUMBER_EDGE_ELEMENTS[i]
                if np.all(Ni == 0):
                    pass
                else:
                    A = Ni[Ni!=0]
                    LEN = len(A)
                    LN = np.arange(start, start+LEN)
                    NUMBER_EDGE_ELEMENTS[i][Ni!=0] = LN
                    start += LEN

        else:
            NUMBER_EDGE_ELEMENTS = None

        NUMBER_EDGE_ELEMENTS = cOmm.bcast(NUMBER_EDGE_ELEMENTS, root=mAster_rank)

        eGM = dict()

        for me in mesh.elements:
            e_MAP = E_MAP[me]
            GV = list()

            for mee in e_MAP:
                GV.append(NUMBER_EDGE_ELEMENTS[mee])

            GV = np.array([_ for _ in chain(*GV)], dtype=int)

            eGM[me] = Gathering_Vector(me, GV)

        eGM = Gathering_Matrix(eGM, mesh_type='_3dCSCG')

        return eGM


    def ___PRIVATE_overcoming_hybrid_singularity_INTERNAL_NODE___(
        self, T, D, C, adt0, e0, Dirichlet_boundaries):
        """"""
        mesh = self._sf_.mesh
        px, py, pz = e0.space.p

        GMt0 = adt0.prime.numbering.gathering
        GMe0 = e0.numbering.gathering

        GMt0_TEW = adt0.prime.numbering.trace_element_wise
        GMe0_EEW = e0.numbering.edge_element_wise

        DICT0= {'WB': 0, 'EB': py, 'WF': -py-1, 'EF': -1,
                'NB': 0, 'SB': px, 'NF': -px-1, 'SF': -1,
                'NW': 0, 'SW': px, 'NE': -px-1, 'SE': -1}

        DICT1= {'W': 0, 'E': -1, 'B': 0, 'F': -1, 'N': 0, 'S':-1}

        DICT3 = {'EF': 'E', 'WF': 'F', 'WB': 'W', 'EB': 'B'}

        sf_local_numbering = self._sf_.numbering.local[0]

        newC = dict()
        newT = dict()

        # if rAnk == mAster_rank:
        #     pbar = tqdm(total=mesh.node.elements.GLOBAL_num, desc='---Parsing-Hybridization-0')

        SKIP_e0_dofs = list() # {e0_dof_number: (edge_element, location), }
        for i in range(mesh.node.elements.GLOBAL_num):

            SOS = mesh.node.elements.do.find.hybrid_singularity_overcoming_setting(i)

            #-------------- INTERNAL NODE SOS ------------------------------------------------------
            if SOS.__class__ == _3dCSCG_InternalNodeSOS:

                #-- We find all the participates: S_e0_dof, N_e0_dof, involved_t0_dofs, involved_e0_dofs --
                trace = SOS.trace

                involved_t0_dofs = dict()
                for _ in trace:
                    TE, CORNER = _
                    if TE in GMt0_TEW:
                        involved_t0_dofs[GMt0_TEW[TE][DICT0[CORNER]]] = _

                edge = SOS.edge
                involved_e0_dofs = dict()
                for _ in edge:
                    EE, END = _
                    if EE in GMe0_EEW:
                        involved_e0_dofs[GMe0_EEW[EE][DICT1[END]]] = _

                if SOS.S_edge in GMe0_EEW:
                    S_e0_dof = GMe0_EEW[SOS.S_edge][0]
                else:
                    S_e0_dof = -1

                if SOS.N_edge in GMe0_EEW:
                    N_e0_dof = GMe0_EEW[SOS.N_edge][-1]
                else:
                    N_e0_dof = -1

                involved_t0_dofs = cOmm.allgather(involved_t0_dofs)
                involved_e0_dofs = cOmm.allgather(involved_e0_dofs)
                S_e0_dof = cOmm.allgather(S_e0_dof)
                N_e0_dof = cOmm.allgather(N_e0_dof)

                ___ = dict()
                for _ in involved_t0_dofs:
                    ___.update(_)
                involved_t0_dofs = ___

                ___ = dict()
                for _ in involved_e0_dofs:
                    ___.update(_)
                involved_e0_dofs = ___

                assert len(involved_t0_dofs) == 4 and len(involved_e0_dofs) == 4, \
                    f"we should have found 4 N-S surface trace dofs and edge dofs."

                for _ in S_e0_dof:
                    if _ != -1:
                        S_e0_dof = _
                        break
                for _ in N_e0_dof:
                    if _ != -1:
                        N_e0_dof = _
                        break

                #- we find the local mesh-elements and indices of all participates ------
                t_dof_mesh_elements, t_dof_local_indices = \
                    GMt0.do.find.elements_and_local_indices_of_dofs(involved_t0_dofs.keys())

                e_dof_mesh_elements, e_dof_local_indices = \
                    GMe0.do.find.elements_and_local_indices_of_dofs(involved_e0_dofs.keys())

                _ = GMe0.do.find.elements_and_local_indices_of_dof(S_e0_dof)
                if _ is not None:
                    Se_mesh_elements, Se_local_indices = _
                else:
                    Se_mesh_elements, Se_local_indices = list(), list()


                _ = GMe0.do.find.elements_and_local_indices_of_dof(N_e0_dof)
                if _ is not None:
                    Ne_mesh_elements, Ne_local_indices = _
                else:
                    Ne_mesh_elements, Ne_local_indices = list(), list()


                # --- clear the existing connection on e0_dofs ------------------------
                for _ in e_dof_mesh_elements:
                    MES = e_dof_mesh_elements[_]
                    IDS = e_dof_local_indices[_]

                    for j, me in enumerate(MES):
                        if me not in newC: newC[me] = C[me].copy().tolil()
                        ind = IDS[j]
                        newC[me][:, ind] = 0 # clear all connections for 0-edge dofs on N-S surface

                #--- Find a first pair ----------------------------------------
                tdf0 = list(involved_t0_dofs.keys())[0]
                position = involved_t0_dofs[tdf0][1]
                cop = DICT3[position]
                for _ in involved_e0_dofs:
                    if involved_e0_dofs[_][1] == cop:
                        break
                edf0 = _
                #>>>>>>>>>>>> connect edf0 to S_e0_dof by tdf0 ----------------------

                tdf_mesh_elements = t_dof_mesh_elements[tdf0]
                tdf_local_indices = t_dof_local_indices[tdf0]

                edf_mesh_elements = e_dof_mesh_elements[edf0]
                edf_local_indices = e_dof_local_indices[edf0]

                for j, me in enumerate(tdf_mesh_elements):
                    if me not in newC: newC[me] = C[me].copy().tolil()
                    if me not in newT: newT[me] = T[me].copy().tolil()

                    ind = tdf_local_indices[j]

                    if me in Se_mesh_elements:
                        # in this mesh-element, we set +1 connection to South-edge-top-dof
                        newT[me][ind, :] = 0
                        assert newC[me][ind].nnz == 0
                        newC[me][ind, Se_local_indices[Se_mesh_elements.index(me)]] = 1

                    else:
                        # in this mesh-element, we set -1 connection to e0_dof_0
                        newT[me][ind, :] = 0
                        assert newC[me][ind].nnz == 0
                        newC[me][ind, edf_local_indices[edf_mesh_elements.index(me)]] = -1

                del involved_t0_dofs[tdf0] # tdf0 will not be used anymore

                #---- Find a second pair ------------------------------------------------------

                for _ in involved_t0_dofs:
                    pos = involved_t0_dofs[_][1]
                    if cop in pos:
                        break
                tdf1 = _

                #>>>>>>>>>>>> connect edf0 to N_e0_dof by tdf1 -----------------------
                tdf_mesh_elements = t_dof_mesh_elements[tdf1]
                tdf_local_indices = t_dof_local_indices[tdf1]

                for j, me in enumerate(tdf_mesh_elements):
                    if me not in newC: newC[me] = C[me].copy().tolil()
                    if me not in newT: newT[me] = T[me].copy().tolil()

                    ind = tdf_local_indices[j]

                    if me in Ne_mesh_elements:
                        # in this mesh-element, we set +1 connection to South-edge-top-dof
                        newT[me][ind, :] = 0
                        assert newC[me][ind].nnz == 0
                        newC[me][ind, Ne_local_indices[Ne_mesh_elements.index(me)]] = 1

                    else:
                        # in this mesh-element, we set -1 connection to e0_dof_1
                        newT[me][ind, :] = 0
                        assert newC[me][ind].nnz == 0
                        newC[me][ind, edf_local_indices[edf_mesh_elements.index(me)]] = -1

                del involved_t0_dofs[tdf1] # tdf1 will not be used anymore
                del involved_e0_dofs[edf0]

                #--- Find a third pair ------------------------------------------------------
                tdf2 = list(involved_t0_dofs.keys())[0]
                position = involved_t0_dofs[tdf2][1]
                cop = DICT3[position]
                for _ in involved_e0_dofs:
                    if involved_e0_dofs[_][1] == cop:
                        break
                edf1 = _

                #>>>>>>>>>>>> connect edf1 to S_e0_dof by tdf2 ----------------------

                tdf_mesh_elements = t_dof_mesh_elements[tdf2]
                tdf_local_indices = t_dof_local_indices[tdf2]

                edf_mesh_elements = e_dof_mesh_elements[edf1]
                edf_local_indices = e_dof_local_indices[edf1]

                for j, me in enumerate(tdf_mesh_elements):
                    if me not in newC: newC[me] = C[me].copy().tolil()
                    if me not in newT: newT[me] = T[me].copy().tolil()

                    ind = tdf_local_indices[j]

                    if me in Se_mesh_elements:
                        # in this mesh-element, we set +1 connection to South-edge-top-dof
                        newT[me][ind, :] = 0
                        assert newC[me][ind].nnz == 0
                        newC[me][ind, Se_local_indices[Se_mesh_elements.index(me)]] = 1

                    else:
                        # in this mesh-element, we set -1 connection to e0_dof_0
                        newT[me][ind, :] = 0
                        assert newC[me][ind].nnz == 0
                        newC[me][ind, edf_local_indices[edf_mesh_elements.index(me)]] = -1

                del involved_t0_dofs[tdf2] # tdf2 will not be used anymore
                del involved_e0_dofs[edf1]

                #---------- The last pair --------------------
                assert len(involved_t0_dofs) == 1
                assert len(involved_e0_dofs) == 2
                for tdf3 in involved_t0_dofs:
                    pass
                position = involved_t0_dofs[tdf3][1]
                cop = DICT3[position]
                for _ in involved_e0_dofs:
                    if involved_e0_dofs[_][1] == cop:
                        break
                edf2 = _

                tdf_mesh_elements = t_dof_mesh_elements[tdf3]
                tdf_local_indices = t_dof_local_indices[tdf3]

                edf_mesh_elements = e_dof_mesh_elements[edf2]
                edf_local_indices = e_dof_local_indices[edf2]

                for j, me in enumerate(tdf_mesh_elements):
                    if me not in newC: newC[me] = C[me].copy().tolil()
                    if me not in newT: newT[me] = T[me].copy().tolil()

                    ind = tdf_local_indices[j]

                    if me in Se_mesh_elements:
                        # in this mesh-element, we set +1 connection to South-edge-top-dof
                        newT[me][ind, :] = 0
                        assert newC[me][ind].nnz == 0
                        newC[me][ind, Se_local_indices[Se_mesh_elements.index(me)]] = 1

                    else:
                        # in this mesh-element, we set -1 connection to e0_dof_0
                        newT[me][ind, :] = 0
                        assert newC[me][ind].nnz == 0
                        newC[me][ind, edf_local_indices[edf_mesh_elements.index(me)]] = -1

                del involved_e0_dofs[edf2]
                #==============================================================================

                for _ in involved_e0_dofs:
                    pass

                if involved_e0_dofs[_][0] in mesh.edge.elements:
                    SKIP_e0_dofs.append(involved_e0_dofs[_])

            #-----------
            elif SOS.__class__ == _3dCSCG_CornerNodeSOS:
                # This `SOS` must only involve one core (all involved mesh-, trace, edge-elements).
                if SOS.mesh is not None:

                    num_trace_elements_on_Dirichlet_boundary = 0

                    for mb in SOS.mesh_boundaries:
                        if mb in Dirichlet_boundaries:
                            num_trace_elements_on_Dirichlet_boundary += 1

                    #------------------- Enough Dirichlet_boundaries found at this corner-node -----
                    if num_trace_elements_on_Dirichlet_boundary >= 2:
                    # this corner-node-element will introduce no singularity if the BC is properly imposed.
                        SKIP_e0_dofs.extend(SOS.edge)

                        tf_DOFs = dict()
                        for ___ in SOS.trace:
                            te, corner = ___[:2]
                            tf_DOFs[GMt0_TEW[te][DICT0[corner]]] = ___

                        corner_mesh_element, corner_name = SOS.mesh

                        s0, s1, s2 = corner_name
                        i0 = 0 if s0 == 'N' else -1
                        i1 = 0 if s1 == 'W' else -1
                        i2 = 0 if s2 == 'B' else -1

                        sf_dof_local_numbering = sf_local_numbering[i0, i1, i2]

                        tf_DOFS_local_indices = list()
                        GV = GMt0[corner_mesh_element]
                        for tdf in tf_DOFs:
                            local_index = GV.index(tdf)
                            tf_DOFS_local_indices.append(local_index)

                        for j, tr in enumerate(SOS.trace):
                            boundary = tr[2]
                            t_local_index = tf_DOFS_local_indices[j]
                            if boundary in Dirichlet_boundaries:
                                if corner_mesh_element in newT:
                                    assert newT[corner_mesh_element][t_local_index].nnz == 0, f"Must be!"
                                else:
                                    assert T[corner_mesh_element][t_local_index].nnz == 0, f"Must be!"
                                assert D[corner_mesh_element][t_local_index].nnz == 1
                                assert D[corner_mesh_element][t_local_index, t_local_index] == 1
                                if corner_mesh_element in newC:
                                    assert newC[corner_mesh_element][t_local_index].nnz == 0
                                else:
                                    assert C[corner_mesh_element][t_local_index].nnz == 0

                            else:

                                assert D[corner_mesh_element][t_local_index].nnz == 0
                                if corner_mesh_element in newT:
                                    assert newT[corner_mesh_element][t_local_index].nnz == 1
                                    assert newT[corner_mesh_element][t_local_index, sf_dof_local_numbering] == 1
                                else:
                                    assert T[corner_mesh_element][t_local_index].nnz == 1
                                    assert T[corner_mesh_element][t_local_index, sf_dof_local_numbering] == 1
                                if corner_mesh_element in newC:
                                    assert newC[corner_mesh_element][t_local_index].nnz == 0
                                else:
                                    assert C[corner_mesh_element][t_local_index].nnz == 0


                    elif num_trace_elements_on_Dirichlet_boundary == 1:
                        # we will have to skip 2-edge dofs
                        ef_DOFs = dict()
                        for ___ in SOS.edge:
                            ee, pos = ___
                            ef_DOFs[GMe0_EEW[ee][DICT1[pos]]] = ___

                        tf_DOFs = dict()
                        for ___ in SOS.trace:
                            te, corner = ___[:2]
                            tf_DOFs[GMt0_TEW[te][DICT0[corner]]] = ___

                        corner_mesh_element, corner_name = SOS.mesh

                        s0, s1, s2 = corner_name
                        i0 = 0 if s0 == 'N' else -1
                        i1 = 0 if s1 == 'W' else -1
                        i2 = 0 if s2 == 'B' else -1
                        sf_dof_local_numbering = sf_local_numbering[i0, i1, i2]

                        if corner_mesh_element not in newC:
                            newC[corner_mesh_element] = C[corner_mesh_element].copy().tolil()
                        if corner_mesh_element not in newT:
                            newT[corner_mesh_element] = T[corner_mesh_element].copy().tolil()

                        ef_DOFS_local_indices = list()
                        GV = GMe0[corner_mesh_element]
                        for edf in ef_DOFs:
                            local_index = GV.index(edf)
                            # clear all existing connections on this edge-dof ----------------
                            newC[corner_mesh_element][:, local_index] = 0
                            ef_DOFS_local_indices.append(local_index)

                        tf_DOFS_local_indices = list()
                        GV = GMt0[corner_mesh_element]
                        for tdf in tf_DOFs:
                            local_index = GV.index(tdf)
                            tf_DOFS_local_indices.append(local_index)

                        NeuB = True
                        for j, tr in enumerate(SOS.trace):
                            boundary = tr[2]
                            t_local_index = tf_DOFS_local_indices[j]
                            if boundary in Dirichlet_boundaries:

                                assert newT[corner_mesh_element][t_local_index].nnz == 0, f"Must be!"
                                assert D[corner_mesh_element][t_local_index, t_local_index] == 1
                                assert D[corner_mesh_element][t_local_index].nnz == 1
                                assert newC[corner_mesh_element][t_local_index].nnz == 0

                                SKIP_e0_dofs.append(SOS.edge[j])

                            else: # on Neumann boundary

                                if NeuB:
                                    assert D[corner_mesh_element][t_local_index].nnz == 0
                                    assert newC[corner_mesh_element][t_local_index].nnz == 0

                                    e_local_index = ef_DOFS_local_indices[j]
                                    newC[corner_mesh_element][t_local_index, e_local_index] = 1

                                    NNZ = newT[corner_mesh_element][t_local_index].nnz
                                    if NNZ == 1:
                                        newT[corner_mesh_element][t_local_index, :] = 0
                                    elif NNZ == 0:
                                        pass
                                    else:
                                        raise Exception()

                                    NeuB = False

                                else:
                                    # we let this pair be
                                    assert D[corner_mesh_element][t_local_index].nnz == 0
                                    assert newC[corner_mesh_element][t_local_index].nnz == 0
                                    NNZ = newT[corner_mesh_element][t_local_index].nnz

                                    if NNZ == 1:
                                        assert newT[corner_mesh_element][t_local_index, sf_dof_local_numbering] == 1
                                    elif NNZ == 0:
                                        newT[corner_mesh_element][t_local_index, sf_dof_local_numbering] = 1
                                    else:
                                        raise Exception()

                                    SKIP_e0_dofs.append(SOS.edge[j])


                    elif num_trace_elements_on_Dirichlet_boundary == 0:
                        # we will have to skip 2-edge dofs
                        ef_DOFs = dict()
                        for ___ in SOS.edge:
                            ee, pos = ___
                            ef_DOFs[GMe0_EEW[ee][DICT1[pos]]] = ___

                        tf_DOFs = dict()
                        for ___ in SOS.trace:
                            te, corner = ___[:2]
                            tf_DOFs[GMt0_TEW[te][DICT0[corner]]] = ___

                        corner_mesh_element, corner_name = SOS.mesh

                        s0, s1, s2 = corner_name
                        i0 = 0 if s0 == 'N' else -1
                        i1 = 0 if s1 == 'W' else -1
                        i2 = 0 if s2 == 'B' else -1
                        sf_dof_local_numbering = sf_local_numbering[i0, i1, i2]

                        if corner_mesh_element not in newC:
                            newC[corner_mesh_element] = C[corner_mesh_element].copy().tolil()
                        if corner_mesh_element not in newT:
                            newT[corner_mesh_element] = T[corner_mesh_element].copy().tolil()

                        ef_DOFS_local_indices = list()
                        GV = GMe0[corner_mesh_element]
                        for edf in ef_DOFs:
                            local_index = GV.index(edf)
                            # clear all existing connections on this edge-dof ----------------
                            newC[corner_mesh_element][:, local_index] = 0
                            ef_DOFS_local_indices.append(local_index)

                        tf_DOFS_local_indices = list()
                        GV = GMt0[corner_mesh_element]
                        for tdf in tf_DOFs:
                            local_index = GV.index(tdf)
                            tf_DOFS_local_indices.append(local_index)

                        for j, _ in enumerate(SOS.trace):

                            t_local_index = tf_DOFS_local_indices[j]

                            assert D[corner_mesh_element][t_local_index].nnz == 0
                            assert newC[corner_mesh_element][t_local_index].nnz == 0
                            NNZ = newT[corner_mesh_element][t_local_index].nnz
                            assert NNZ <= 1

                            if j == 0:
                                if NNZ == 1:
                                    assert newT[corner_mesh_element][t_local_index, sf_dof_local_numbering] == 1
                                else:
                                    newT[corner_mesh_element][t_local_index, sf_dof_local_numbering] = 1

                            else:
                                e_local_index = ef_DOFS_local_indices[j]
                                if NNZ == 1:
                                    newT[corner_mesh_element][t_local_index, sf_dof_local_numbering] = 0

                                newC[corner_mesh_element][t_local_index, e_local_index] = 1

                    else:
                        raise Exception()

            # elif SOS.__class__ == _3dCSCG_CornerEdgeNodeSOS:
            #     pass
            # elif SOS.__class__ == _3dCSCG_SurfaceMiddleNodeSOS:
            #     pass
            # else:
            #     raise NotImplementedError(f"Cannot handle SOS: {SOS.__class__.__name__}")

        #     #============================================================================
        #     if rAnk == mAster_rank:
        #         pbar.update(1)
        # #=====================================================================================
        # if rAnk == mAster_rank:
        #     pbar.close()

        for _ in T:
            if _ not in newT:
                newT[_] = T[_]
            else:
                # noinspection PyUnresolvedReferences
                newT[_] = newT[_].tocsr()

        for _ in C:
            if _ not in newC:
                newC[_] = C[_]
            else:
                # noinspection PyUnresolvedReferences
                newC[_] = newC[_].tocsr()

        newT = T.__class__(mesh, newT, cache_key_generator = 'no_cache')
        newC = C.__class__(mesh, newC, cache_key_generator = 'no_cache')

        return newT, newC, SKIP_e0_dofs




if __name__ == '__main__':
    # mpiexec -n 5 python objects\CSCG\_3d\forms\standard\_0s\special\main.py
    from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller, ExactSolutionSelector

    elements = [2, 3, 2]
    mesh = MeshGenerator('crazy_periodic', c=0.1)(elements)
    # mesh = MeshGenerator('crazy', c=0.1)(elements)
    # mesh = MeshGenerator('cuboid_periodic', region_layout=[3,2,2])(elements)
    # mesh = MeshGenerator('cuboid', region_layout=[3,2,2])(elements)

    ES = ExactSolutionSelector(mesh)('Poisson:sincos1')

    Dirichlet_boundaries = [] #
    Neumann_boundaries = []
    # Dirichlet_boundaries = ['Front', 'West', "East", 'Back', ] #
    # Neumann_boundaries = ['North', 'South', ]

    space = SpaceInvoker('polynomials')([2, 1, 3])
    FC = FormCaller(mesh, space)

    f = FC('0-f', is_hybrid=True)
    t = FC('0-adt')
    e = FC('0-e')



    f.TW.BC.body = ES.status.potential
    f.TW.do.push_BC_to_instant(0)
    f.BC.valid_boundaries = Neumann_boundaries

    # t.prime.TW.BC.body = ES.status.velocity.flux
    t.prime.TW.BC.body = ES.status.potential
    t.prime.TW.do.push_BC_to_instant(0)
    t.BC.valid_boundaries = Dirichlet_boundaries

    T, D, C, b, eGM = f.special.hybrid_pairing(t, e)

    # f.TW.BC.body = ES.status.velocity
    # f.TW.do.push_BC_to_instant(0)
    # f.BC.valid_boundaries = Neumann_boundaries
    #
    # t.prime.TW.BC.body = ES.status.velocity.components.T_perp
    # t.prime.TW.do.push_BC_to_instant(0)
    # t.BC.valid_boundaries = Dirichlet_boundaries
    #
    # T, D, C, b, eGM = f.special.hybrid_pairing(t, e)

    # T = t.matrices.trace
    # C = e.matrices.complement
    # T, C = f.special.___PRIVATE_overcoming_hybrid_singularity___(
    #     T, C, t, e, Dirichlet_boundaries=Dirichlet_boundaries)

    # f.dofs.visualize.matplot.connection_around_node_element(7, T, D, C, t, e, checking_mode=False)

    for i in range(mesh.node.elements.GLOBAL_num):
        f.dofs.visualize.matplot.connection_around_node_element(i, T, D, C, t, e, checking_mode=True)

    # for i in range(t1.prime.numbering.gathering.GLOBAL_num_dofs):
    #     f1.dofs.visualize.matplot.connection_through_trace_dof(i, T, C, t1, e1, checking_mode=True)
    #
    # for i in range(e1.numbering.gathering.GLOBAL_num_dofs):
    #     f1.dofs.visualize.matplot.connection_through_around_edge_dof(
    #         i, T, C, t1, e1, checking_mode=True)