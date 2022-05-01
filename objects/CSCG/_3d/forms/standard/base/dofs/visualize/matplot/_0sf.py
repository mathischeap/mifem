import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.main import FrozenOnly
from screws.miscellaneous.matrix_simplify import rsmat
from root.config.main import rAnk, mAster_rank, cOmm, np, MPI

import matplotlib.pyplot as plt
import matplotlib
from scipy.sparse import csr_matrix




class _3dCSCG_S0F_DOFs_Matplot(FrozenOnly):
    """"""
    def __init__(self, dofs):
        """"""
        self._dofs_ = dofs
        self._sf_ = dofs._sf_
        self._mesh_ = dofs._sf_.mesh
        self._freeze_self_()

    def __call__(self, *args, **kwargs):
        """"""
        raise NotImplementedError()

    def dofs_around_node_element(self, i, zoom=0.9):
        """We plot the dofs around the node-element #i.

        Parameters
        ----------
        i :
        zoom :

        Returns
        -------

        """
        assert self._sf_.IS.hybrid, f"Visualize dofs around a node-element is only valid for " \
                                    f"HYBRID 0-standard-form."

        #----------- plot the trace-elements and mesh-elements around node-element #i --------
        if i in self._mesh_.node.elements:
            node = self._mesh_.node.elements[i]
            positions = node.positions
            DATA = dict()
            node_LOCATION = list()
            TraceElement_Data = dict()

            on_mesh_boundary = node.IS.on_mesh_boundary

            TE_coordinates = (
                [-1, 0, 0],
                [ 1, 0, 0],
                [0, -1, 0],
                [0,  1, 0],
                [0, 0, -1],
                [0, 0,  1],
            )

            for pos in positions:
                if pos[0] in '0123456789': # at the corner of a mesh element
                    element_num = int(pos[:-3])
                    if element_num in self._mesh_.node.elements._mesh_.elements:
                        element = self._mesh_.node.elements._mesh_.elements[element_num]
                        data = element.do.generate_element_plot_data(zoom=zoom)
                        DATA[element_num] = (data, pos[-3:])

                        node_LOCATION .append(
                            node.coordinate_transformation.mapping(
                            from_element=element_num,
                                corner=pos[-3:])
                        )

                        TE_MAP = self._mesh_.node.elements._mesh_.trace.elements.map[element_num]

                        TraceElement_Data[element_num] = list()

                        IND = ['NSWEBF'.index(_) for _ in pos[-3:]]

                        for _, t in enumerate(TE_MAP):
                            if _ in IND:
                                coo = element.coordinate_transformation.mapping(*TE_coordinates[_])
                                TraceElement_Data[element_num].append([t, coo])

        else:
            positions = None
            DATA = dict()
            node_LOCATION = None
            on_mesh_boundary = None
            TraceElement_Data = dict()

        positions = cOmm.allgather(positions)
        on_mesh_boundary = cOmm.allgather(on_mesh_boundary)
        for _ in positions:
            if _ is not None:
                positions = _
                break

        DATA = cOmm.gather(DATA, root=mAster_rank)
        TraceElement_Data = cOmm.gather(TraceElement_Data, root=mAster_rank)
        node_LOCATION = cOmm.gather(node_LOCATION, root=mAster_rank)
        on_mesh_boundary = cOmm.reduce(on_mesh_boundary, root=mAster_rank, op=MPI.LOR)

        nL = list()
        if rAnk == mAster_rank:
            ___ = dict()
            for data in DATA:
                ___.update(data)
            DATA = ___

            ___ = dict()
            for data in TraceElement_Data:
                ___.update(data)
            TraceElement_Data = ___

            for __ in node_LOCATION:
                if __ is not None:
                    for _ in __:
                        nL.append(_)

            fig = plt.figure(figsize=(12, 8))
            ax = fig.add_subplot(111, projection='3d')
            # make the panes transparent
            ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
            ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
            ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
            # make the grid lines transparent
            ax.xaxis._axinfo["grid"]['color'] =  (1,1,1,0)
            ax.yaxis._axinfo["grid"]['color'] =  (1,1,1,0)
            ax.zaxis._axinfo["grid"]['color'] =  (1,1,1,0)

            x_lim, y_lim, z_lim = [list() for _ in range(3)]

            for nl in nL:
                ax.scatter(*nl, marker='o', facecolor='w', color='b')

            for e in DATA:
                data, pos = DATA[e]
                X = data[:,0,:]
                Y = data[:,1,:]
                Z = data[:,2,:]
                x_lim.append(np.max(X))
                x_lim.append(np.min(X))
                y_lim.append(np.max(Y))
                y_lim.append(np.min(Y))
                z_lim.append(np.max(Z))
                z_lim.append(np.min(Z))

                for line in data:
                    ax.plot(*line, color='gray', linewidth=0.75)

                ax.text(np.mean(X), np.mean(Y), np.mean(Z), str(e),
                            bbox=dict(boxstyle="square",
                                       ec=(0., 1, 1),
                                       fc=(0., 1, 1),
                                       ),
                        ha='center', va='center', ma='center')

                Ts_COOs = TraceElement_Data[e]
                for t_coo in Ts_COOs:
                    t, coo = t_coo
                    cx, cy, cz = coo
                    cx = cx[0]
                    cy = cy[0]
                    cz = cz[0]
                    ax.text(cx, cy, cz, str(t), color='g',
                            bbox=dict(boxstyle="square",
                                       ec=(1., 0.5, 0.5, 0.5),
                                       fc=(1., 0.8, 0.8, 0.25),
                                       ),
                            ha='center', va='center', ma='center')

            ax.set_box_aspect((np.ptp(x_lim), np.ptp(y_lim), np.ptp(z_lim)))

            ax.set_xlabel(r'$x$', fontsize=12)
            ax.set_ylabel(r'$y$', fontsize=12)
            ax.set_zlabel(r'$z$', fontsize=12)
            if on_mesh_boundary:
                plt.title(f'Around node-element #{i} on mesh boundary.')
            else:
                plt.title(f'Around node-element #{i}')

        else:
            pass

        #----------------- we now to fond all form dofs around node-element------------------------------
        sf_DOFS = list()
        for pos in positions:
            if pos[0].isnumeric():  # we found a mesh-element position
                mesh_element = int(pos[:-3])
                corner = pos[-3:]
                dof = self._dofs_.do.find.dof_at_corner_of_mesh_element(mesh_element, corner)
                if dof is not None:
                    sf_DOFS.append(dof)

        sf_DOFS = cOmm.allgather(sf_DOFS)

        __ = set()
        for _ in sf_DOFS:
            __.update(_)
        sf_DOFS = list(__)
        sf_DOFS.sort()

        #------------------ find the coordinates of all sf dofs ---------------------------
        if rAnk == mAster_rank:
            ALL_sf_dofs_COO = dict()
        for sfd in sf_DOFS:
            dof = self._dofs_[sfd]
            COO = dof.do.generate_plot_data_from_element(None, zoom=zoom)['DOF_COO']

            COO = cOmm.gather(COO, root=mAster_rank)

            if rAnk == mAster_rank:
                ___ = dict()
                for __ in COO:
                    ___.update(__)
                COO = ___

                ALL_sf_dofs_COO[sfd] = COO

        #-------------------------------------------------------------------------------------------

        if rAnk == mAster_rank:
            for sfd in ALL_sf_dofs_COO:
                all_positions = ALL_sf_dofs_COO[sfd]
                for me in all_positions:
                    x, y, z = all_positions[me]
                    ax.scatter(x, y, z, marker='s')
            plt.show()
            plt.close(fig)

    def connection_around_node_element(
        self, i, T, D, C, tf, ef,
        checking_mode=False, zoom=0.85, density=5,
        saveto=None, title=None):
        """

        Parameters
        ----------
        i
        T :
            Trace matrix of the `tf`, a 0-(AD)trace-form.
        D :
        C :
            The complement matrix of the ef, a 0-edge-form.
        tf :
            The 0-(AD)trace-form.
        ef :
            The 0-edge-form.
        checking_mode : bool, optional
            {default: False}
        zoom : float, optional
        density : int, optional
        saveto : bool, optional
        title : {None, str}, optional

        Returns
        -------

        """
        assert self._sf_.IS.hybrid, f"the form is not hybrid, this method makes no sense."

        assert T.__class__.__name__ == 'EWC_SparseMatrix', f"T must be a trace matrix."
        assert C.__class__.__name__ == 'EWC_SparseMatrix', f"C must be a complement matrix."
        mesh = self._sf_.mesh

        if hasattr(tf, 'prime'): tf = tf.prime
        assert tf.__class__.__name__ == '_3dCSCG_0Trace'
        assert ef.__class__.__name__ == '_3dCSCG_0Edge'

        assert self._sf_.mesh == ef.mesh
        assert self._sf_.mesh == tf.mesh
        assert self._sf_.space == ef.space
        assert self._sf_.space == tf.space

        sGM = self._sf_.numbering.gathering
        tGM = tf.numbering.gathering
        eGM = ef.numbering.gathering
        # ST = tf.matrices.selective.T

        sf = self._sf_

        # the coordinate of this node-element -----------------------------
        if i in mesh.node.elements:
            node_element = mesh.node.elements[i]

            if node_element.IS.on_periodic_boundary:

                N_xyz = list()
                positions = node_element.positions
                for pos in positions:
                    if pos[0].isnumeric():
                        mesh_element = int(pos[:-3])

                        if mesh_element in mesh.elements:

                            N_xyz .append(
                                node_element.coordinate_transformation.mapping(
                                    from_element=mesh_element,
                                    corner = pos[-3:])
                            )

            else:
                N_xyz = [node_element.coordinate_transformation.mapping(),]
        else:
            N_xyz = None

        N_xyz = cOmm.gather(N_xyz, root=mAster_rank)

        #--------- bcast positions, on_mesh_boundary to all cores ----------------------------------------------
        if i in mesh.node.elements:
            node_element = mesh.node.elements[i]
            positions = node_element.positions
            on_mesh_boundary = node_element.IS.on_mesh_boundary

        else:
            positions = None
            on_mesh_boundary = False

        positions = cOmm.allgather(positions)
        on_mesh_boundary = cOmm.allgather(on_mesh_boundary)

        for _ in positions:
            if _ is not None:
                positions = _
                break
        on_mesh_boundary = any(on_mesh_boundary)


        # plot data of the mesh-elements-around the node-element --------------------------
        MEF = dict()
        for pos in positions:
            if pos[0].isnumeric(): # we found a mesh-element position
                mesh_element = int(pos[:-3])
                if mesh_element in mesh.elements: # this is a local mesh element
                    data = mesh.elements[mesh_element].do.generate_element_plot_data(
                        zoom=zoom, density=density)
                    MEF[mesh_element] = data
            else:
                pass
        MEF = cOmm.gather(MEF, root=mAster_rank)

        # find all 0-standard-form dofs at this node element ---------------------------
        sf_DOFS = list()
        for pos in positions:
            if pos[0].isnumeric(): # we found a mesh-element position
                mesh_element = int(pos[:-3])
                corner = pos[-3:]
                dof = self._dofs_.do.find.dof_at_corner_of_mesh_element(mesh_element, corner)
                if dof is not None:
                    sf_DOFS.append(dof)

        sf_DOFS = cOmm.allgather(sf_DOFS)

        __ = set()
        for _ in sf_DOFS:
            __.update(_)
        sf_DOFS = list(__)
        sf_DOFS.sort()

        #--------- find all 0-(AD)trace-form dofs at this node element ------------------------
        tf_DOFS = list()
        for pos in positions:
            if pos[0].isnumeric(): # we found a mesh-element position
                mesh_element = int(pos[:-3])
                corner = pos[-3:]
                dof = tf.dofs.do.find.dof_at_corner_of_mesh_element(mesh_element, corner)
                if dof is not None:
                    tf_DOFS.extend(dof)

        tf_DOFS = cOmm.allgather(tf_DOFS)

        __ = set()
        for _ in tf_DOFS:
            __.update(_)
        tf_DOFS = list(__)
        tf_DOFS.sort()

        #--------- find all 0-edge-form dofs at this node element -----------------------------
        ef_DOFS = list()
        for pos in positions:
            if pos[0].isnumeric(): # we found a mesh-element position
                mesh_element = int(pos[:-3])
                corner = pos[-3:]
                dof = ef.dofs.do.find.dof_at_corner_of_mesh_element(mesh_element, corner)
                if dof is not None:
                    ef_DOFS.extend(dof)

        ef_DOFS = cOmm.allgather(ef_DOFS)

        __ = set()
        for _ in ef_DOFS:
            __.update(_)
        ef_DOFS = list(__)
        ef_DOFS.sort()

        #------------------ find the coordinates of all sf dofs ---------------------------
        if rAnk == mAster_rank:
            ALL_sf_dofs_COO = dict()
        for sfd in sf_DOFS:
            dof = sf.dofs[sfd]
            COO = dof.do.generate_plot_data_from_element(None, zoom=zoom)['DOF_COO']

            COO = cOmm.gather(COO, root=mAster_rank)

            if rAnk == mAster_rank:
                ___ = dict()
                for __ in COO:
                    ___.update(__)
                COO = ___

                ALL_sf_dofs_COO[sfd] = COO

        #============ STARTING CHECKING ======================================================
        sf_ef_dof_dict = dict()
        I = 0
        for dof in sf_DOFS:
            sf_ef_dof_dict['s-' + str(dof)] = I
            I += 1
        for dof in ef_DOFS:
            sf_ef_dof_dict['e-' + str(dof)] = I
            I += 1

        tf_dof_dict = dict()
        I = 0
        for dof in tf_DOFS:
            tf_dof_dict['t-' + str(dof)] = I
            I += 1

        if not on_mesh_boundary: assert len(tf_dof_dict) == 12 and len(sf_ef_dof_dict) == 14

        MATRIX = np.zeros((len(tf_dof_dict), len(sf_ef_dof_dict)))

        for tdf in tf_DOFS:
            _ = tGM.do.find.elements_and_local_indices_of_dof(tdf)
            if _ is not None:

                tf_indicator = 't-' + str(tdf)

                for ele, ind in zip(*_):
                    T_row = T[ele][ind].tocsr()
                    C_row = C[ele][ind].tocsr()

                    for tri, trd in zip(T_row.indices, T_row.data):
                        # involved sf dof: ele, tri
                        involved_sf_dof = sGM[ele][tri]
                        assert involved_sf_dof in sf_DOFS, "must be the case!"
                        sf_indicator = 's-' + str(involved_sf_dof)

                        MATRIX[tf_dof_dict[tf_indicator], sf_ef_dof_dict[sf_indicator]] = trd

                    for cri, crd in zip(C_row.indices, C_row.data):
                        # involved sf dof: ele, tri
                        involved_ef_dof = eGM[ele][cri]
                        assert involved_ef_dof in ef_DOFS, "must be the case!"
                        ef_indicator = 'e-' + str(involved_ef_dof)

                        MATRIX[tf_dof_dict[tf_indicator], sf_ef_dof_dict[ef_indicator]] = crd

                    if T_row.nnz == 1 and C_row.nnz == 0:
                        assert D[ele][ind].nnz == 0
                    elif T_row.nnz == 0 and C_row.nnz == 0:
                        assert D[ele][ind].nnz == 1
                    elif T_row.nnz == 0 and C_row.nnz == 1:
                        assert D[ele][ind].nnz == 0
                    else:
                        raise NotImplementedError()

            if not on_mesh_boundary:
                if _ is not None:
                    for ele, ind in zip(*_):
                        assert D[ele][ind].nnz == 0

        MATRIX = cOmm.gather(MATRIX, root=mAster_rank)
        if rAnk == mAster_rank:

            EMPTY_tf_dofs = list()

            MATRIX = np.sum(MATRIX, axis=0)

            if not on_mesh_boundary: # an internal node-element

                I = 0
                for j in range(14):
                    Mj = MATRIX[:, j]
                    if np.all(Mj==0):
                        I += 1

                assert I == 1, f"For an internal node-element, we left one 0-edge form dof to be" \
                               f"unused!"

                for Mi in MATRIX:
                    Mi = csr_matrix(Mi)
                    assert Mi.nnz == 2
                    i0, i1 = Mi.data
                    if i0 == 1:
                        assert i1 == -1
                    elif i0 == -1:
                        assert i1 == 1
                    else:
                        raise Exception()

                MATRIX = rsmat(MATRIX)
                MATRIX = MATRIX[-1,:]
                if np.all(MATRIX == 0):
                    raise Exception(f"Singularity find around internal node-element#{i}.")
            else:

                NEW_MATRIX = list()

                for j, Mi in enumerate(MATRIX):

                    if csr_matrix(Mi).nnz == 0:
                        for tf_indicator in tf_dof_dict:
                            if tf_dof_dict[tf_indicator] == j:
                                break
                        tf_dof = int(tf_indicator[2:])
                        EMPTY_tf_dofs.append(tf_dof)

                    else:
                        NEW_MATRIX.append(Mi)

                NEW_MATRIX = np.array(NEW_MATRIX)

                if NEW_MATRIX.shape[0] == 0:
                    pass
                else:
                    NEW_MATRIX = rsmat(NEW_MATRIX)
                    NEW_MATRIX = NEW_MATRIX[-1,:]
                    if np.all(NEW_MATRIX == 0):
                        raise Exception(f"Singularity find around boundary node-element#{i}.")

        else:
            EMPTY_tf_dofs = None

        EMPTY_tf_dofs = cOmm.bcast(EMPTY_tf_dofs, root=mAster_rank)

        for tdf in EMPTY_tf_dofs:
            _ = tGM.do.find.elements_and_local_indices_of_dof(tdf)
            if _ is not None:
                for ele, ind in zip(*_):
                    assert D[ele][ind, ind] == 1

        #============ BELOW::: DO THE PLOT ====================================================
        if rAnk != mAster_rank: return
        if checking_mode: return

        if saveto is not None: matplotlib.use('Agg')
        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(111, projection='3d')
        # make the panes transparent
        ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        # make the grid lines transparent
        ax.xaxis._axinfo["grid"]['color'] = (1, 1, 1, 0)
        ax.yaxis._axinfo["grid"]['color'] = (1, 1, 1, 0)
        ax.zaxis._axinfo["grid"]['color'] = (1, 1, 1, 0)
        ax.tick_params(labelsize=12)
        ax.set_xlabel(r'$x$', fontsize=15)
        ax.set_ylabel(r'$y$', fontsize=15)
        ax.set_zlabel(r'$z$', fontsize=15)

        # ------- plot all mesh elements around the node element -----------------------------
        for ___ in MEF:
            for me in ___:
                efd = ___[me]
                for line in efd:
                    ax.plot(*line, color='gray', linewidth=0.75)

        #--------- scatter all sf dofs -------------------------------------------------------
        for sfd in ALL_sf_dofs_COO:
            all_positions = ALL_sf_dofs_COO[sfd]
            for me in all_positions:
                x, y, z = all_positions[me]
                ax.scatter(x, y, z, marker='s')

        #---------- SCATTER THE NODE-ELEMENT ------------------------------------------------
        for XYZ in N_xyz:
            if XYZ is not None:
                for xyz in XYZ:
                    x, y, z = xyz
                    x, y, z = x[0], y[0], z[0]
                    ax.scatter(x, y, z, color='b', facecolor='w')

        # --------- title ------------------------------------------------------------------
        if title is None:
            plt.title(f"Connection-around-1-node-element#{i}.")
        else:
            plt.title(title)

        if saveto is None:
            plt.show()
        else:
            plt.savefig(saveto, bbox_inches='tight')

        plt.close()






if __name__ == '__main__':
    # mpiexec -n 5 python objects\CSCG\_3d\forms\standard\base\dofs\visualize\matplot\_0sf.py
    from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller

    elements = [8,8,8]
    # mesh = MeshGenerator('crazy_periodic', c=0.1)(elements)
    mesh = MeshGenerator('crazy', c=0.1)(elements)
    # Dirichlet_boundaries = []
    Dirichlet_boundaries = ['North', 'South', 'West', 'East', 'Back', 'Front',] #
    # mesh = MeshGenerator('bridge_arch_cracked')(elements)
    space = SpaceInvoker('polynomials')([2, 2, 2])
    FC = FormCaller(mesh, space)

    f = FC('0-f', is_hybrid=True)
    t = FC('0-adt')
    e = FC('0-e')

    T = t.matrices.trace
    C = e.matrices.complement
    T, C = f.special.___PRIVATE_overcoming_hybrid_singularity___(
        T, C, t, e, Dirichlet_boundaries=Dirichlet_boundaries)

    f.dofs.visualize.matplot.connection_around_node_element(7, T, C, t, e, checking_mode=False)

    # for i in range(mesh.edge.elements.GLOBAL_num):
    #     f.dofs.visualize.matplot.dofs_around_node_element(i)