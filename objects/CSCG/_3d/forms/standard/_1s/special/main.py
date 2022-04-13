
import sys
if './' not in sys.path: sys.path.append('./')

from screws.freeze.main import FrozenOnly
from tools.linear_algebra.elementwise_cache.objects.sparse_matrix.main import EWC_SparseMatrix
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
        (self X 1-form, 1-form). To first cross product with a 1-form then do a inner product with
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



    def overcoming_hybrid_singularity_with_BC(self, T, C, Dirichlet_boundaries):
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
        full_C = self.overcoming_hybrid_singularity(T, C, Dirichlet_boundaries=None)[1]
        new_CT_full = full_C.T


        new_T = dict() # the new Trace matrix.
        new_C = dict() # the new Trace matrix.

        DB_Block = dict()

        return new_T, DB_Block, new_C, new_CT_full






    def overcoming_hybrid_singularity(self, T, C, Dirichlet_boundaries=None):
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

        T_shape = T.shape
        C_shape = C.shape
        k = self._sf_.k

        num_T_basis = getattr(self._sf_.space.num_basis, f'_3dCSCG_{k}Trace')[0]
        num_S_basis = self._sf_.num.basis
        num_E_basis = getattr(self._sf_.space.num_basis, f'_3dCSCG_{k}Edge')[0]
        mesh = self._sf_.mesh
        num_elements = len(mesh.elements)

        assert T_shape == (num_elements, num_T_basis, num_S_basis), f"T shape does not match."
        assert C_shape == (num_elements, num_T_basis, num_E_basis), f"C shape does not match."

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
                pass

            else:

                SOS = mesh.edge.elements.do.find.hybrid_singularity_overcoming_setting(i)

                if SOS is None: # this core has no business with this SOS
                    pass
                else:
                    i, replacing, through = SOS.i, SOS.replacing, SOS.through
                    mesh_element, corner_edge = replacing
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

        return nT, nC







if __name__ == '__main__':
    # mpiexec -n 5 python objects\CSCG\_3d\forms\standard\_1s\special\main.py
    from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller

    elements = [2,2,2]
    # mesh = MeshGenerator('crazy_periodic', c=0.1)(elements)
    mesh = MeshGenerator('crazy', c=0.1)(elements)

    # Dirichlet_boundaries = ['North', 'South', 'West', 'East', 'Back', 'Front',] #
    # mesh = MeshGenerator('bridge_arch_cracked')(elements)
    space = SpaceInvoker('polynomials')([2, 2, 2])
    FC = FormCaller(mesh, space)

    f1 = FC('1-f', is_hybrid=True)
    t1 = FC('1-adt')
    e1 = FC('1-e')

    T = t1.matrices.trace
    C = e1.matrices.complement

    T, C = f1.special.overcoming_hybrid_singularity(T, C)

    # f1.dofs.visualize.matplot.connection_through_trace_dof(55, T, C, t1, e1, checking_mode=True)


    #
    for i in range(t1.prime.numbering.gathering.GLOBAL_num_dofs):
        f1.dofs.visualize.matplot.connection_through_trace_dof(i, T, C, t1, e1, checking_mode=False)
