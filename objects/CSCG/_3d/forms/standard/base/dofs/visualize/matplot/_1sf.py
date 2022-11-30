# -*- coding: utf-8 -*-

import sys
if './' not in sys.path: sys.path.append('./')
import matplotlib.pyplot as plt
import matplotlib

from components.freeze.main import FrozenOnly
from root.config.main import COMM, RANK, MASTER_RANK
from components.generators.counter import Counter
import numpy as np
from components.miscellaneous.matrix_simplify import rsmat
from components.warnings.hybrid_singularity import HybridSingularityWarning
import warnings



class _3dCSCG_S1F_DOFs_Matplot(FrozenOnly):
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

    def connection_through_trace_dof(self,
        i, T, C, tf, ef,
        density=5, linewidth=0.75, zoom=0.9, title=None, checking_mode=False,
        saveto=None,
        ):
        """

        Parameters
        ----------
        i : int
            the number of a dof of a 1-(AD)trace form.
        T :
            The trace matrix.
        C :
            The complementary matrix.
        tf :
            The (AD)trace form.
        ef :
            The edge form.
        density
        linewidth
        zoom
        title
        checking_mode : bool, optional
            {default: False} If `checking_mode` is on (True), we only skip the plotting (to
            accelerate the program.)
        saveto

        Returns
        -------

        """
        assert self._sf_.whether.hybrid, f"the form is not hybrid, this method makes no sense."

        assert T.__class__.__name__ == 'EWC_SparseMatrix', f"T must be a trace matrix."
        assert C.__class__.__name__ == 'EWC_SparseMatrix', f"C must be a complement matrix."
        mesh = self._sf_.mesh

        if hasattr(tf, 'prime'): tf = tf.prime
        assert tf.__class__.__name__ == '_3dCSCG_1Trace'
        assert ef.__class__.__name__ == '_3dCSCG_1Edge'

        assert self._sf_.mesh == ef.mesh
        assert self._sf_.mesh == tf.mesh
        assert self._sf_.space == ef.space
        assert self._sf_.space == tf.space

        GMt = tf.numbering.gathering
        GMe = ef.numbering.gathering
        GMs = self._sf_.numbering.gathering

        #---------------- trace element i plotting data ------------------------------------

        TEi_PD = tf.dofs[i].do.generate_plot_data(density=2*density, zoom=zoom)

        TEi_PD = COMM.gather(TEi_PD, root=MASTER_RANK)
        if RANK == MASTER_RANK:
            ___ = list()
            for _ in TEi_PD:
                if _ is not None:
                    ___.append(_)
            TEi_PD = ___

            Trace_Element_Frames = list()
            Trace_dof_xyz = list()
            for _ in TEi_PD:
                Trace_Element_Frames.append(_['Trace-Element-Frame'])
                Trace_dof_xyz.append(_['DOF-xyz'])

        #----------------------------------------------------------------------------------------
        t_elements_indices = GMt.do.find.elements_and_local_indices_of_dof(i)

        if t_elements_indices is not None:

            INVOLVED_sf_local_dofs = list()
            INVOLVED_ef_local_dofs = list()

            t_elements, t_indices = t_elements_indices
            for te, ti in zip(t_elements, t_indices):

                assert te in mesh.elements

                sf_row = T[te][ti]
                if sf_row.__class__.__name__ == 'csc_matrix':
                    sf_row = sf_row.tocsr()
                data, indices = sf_row.data, sf_row.indices

                for d, I in zip(data, indices):
                    INVOLVED_sf_local_dofs.append((te, I, GMs[te][I], d))


                ef_row = C[te][ti]
                if ef_row.__class__.__name__ == 'csc_matrix':
                    ef_row = ef_row.tocsr()

                data, indices = ef_row.data, ef_row.indices
                for d, I in zip(data, indices):
                    INVOLVED_ef_local_dofs.append((te, I, GMe[te][I], d))

        else:
            INVOLVED_sf_local_dofs = None
            INVOLVED_ef_local_dofs = None

        #---------------------------------------------------------------------------------------
        INVOLVED_sf_local_dofs = COMM.allgather(INVOLVED_sf_local_dofs)
        INVOLVED_ef_local_dofs = COMM.allgather(INVOLVED_ef_local_dofs)

        ___ = list()
        for _ in INVOLVED_sf_local_dofs:
            if _ is not None:
                ___.extend(_)
        INVOLVED_sf_local_dofs = ___
        # [(1, 7, 19, 1), (in this mesh-element, local_index, global numbering, value (+-1)),]

        ___ = list()
        for _ in INVOLVED_ef_local_dofs:
            if _ is not None:
                ___.extend(_)
        INVOLVED_ef_local_dofs = ___

        assert len(INVOLVED_ef_local_dofs) <= 1
        if len(INVOLVED_ef_local_dofs) == 1:
            if len(INVOLVED_sf_local_dofs) == 1:
                v_sf = INVOLVED_sf_local_dofs[0][3]
                v_ef = INVOLVED_ef_local_dofs[0][3]
                assert v_sf + v_ef == 0, "Must form a pair."
            else:
                assert len(INVOLVED_sf_local_dofs) == 0, f"must be this case"

        #---------------------------------------------------------------------------------------
        EFD = list()
        XYZ = list()
        for sf_dof in INVOLVED_sf_local_dofs:
            mesh_element, local_numbering, global_numbering, value = sf_dof
            dof = self._dofs_[global_numbering]

            gpd = dof.do.generate_plot_data_from_element(mesh_element, density=2*density, zoom=zoom)
            efd = gpd['EFD']
            xyz = gpd['DOF_COO']

            if efd is not None: EFD.append(efd)
            if xyz is not None: XYZ.append([xyz, value])

        #---------------------------------------------------------------------------------------

        if len(INVOLVED_ef_local_dofs) == 1:
            mesh_element, local_numbering, global_numbering, value = INVOLVED_ef_local_dofs[0]
            eEFD = ef.special.generate_plot_data_for_dof(
                global_numbering, density=density, zoom=zoom)
        else:
            eEFD = None

        EFD = COMM.gather(EFD, root=MASTER_RANK)
        XYZ = COMM.gather(XYZ, root=MASTER_RANK)

        #----------------------------------------------------------------------------------------
        if RANK == MASTER_RANK:

            ___ = list()
            for _ in XYZ: ___.extend(_)
            XYZ = ___

            v1_v2 = list()
            for xyz_v in XYZ:
                _, v = xyz_v
                v1_v2.append(v)

            if len(v1_v2) == 2: # we find two standard-form-dofs, they must cancel!
                assert np.sum(v1_v2) == 0, 'two 1-form dof must enter with different sign!'

            if checking_mode:
                return

            #------------ plot ----------------------------------------------------------------
            ___ = list()
            for _ in EFD: ___.extend(_)
            EFD = ___

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

            for efd in EFD:
                if 'xLines_x' in efd:
                    X, Y, Z = efd['xLines_x'], efd['xLines_y'], efd['xLines_z']
                    for x, y, z in zip(X, Y, Z): plt.plot(x, y, z, 'lightgray', linewidth=linewidth)
                if 'yLines_x' in efd:
                    X, Y, Z = efd['yLines_x'], efd['yLines_y'], efd['yLines_z']
                    for x, y, z in zip(X, Y, Z): plt.plot(x, y, z, 'lightgray', linewidth=linewidth)
                if 'zLines_x' in efd:
                    X, Y, Z = efd['zLines_x'], efd['zLines_y'], efd['zLines_z']
                    for x, y, z in zip(X, Y, Z): plt.plot(x, y, z, 'lightgray', linewidth=linewidth)

                X, Y, Z = efd['xLines_x_B'], efd['xLines_y_B'], efd['xLines_z_B']
                for x, y, z in zip(X, Y, Z): plt.plot(x, y, z, 'green', linewidth=linewidth)
                X, Y, Z = efd['yLines_x_B'], efd['yLines_y_B'], efd['yLines_z_B']
                for x, y, z in zip(X, Y, Z): plt.plot(x, y, z, 'green', linewidth=linewidth)
                X, Y, Z = efd['zLines_x_B'], efd['zLines_y_B'], efd['zLines_z_B']
                for x, y, z in zip(X, Y, Z): plt.plot(x, y, z, 'green', linewidth=linewidth)

                x, y, z = efd['center']['coordinate']
                x, y, z = x[0], y[0], z[0]
                num = efd['center']['number']
                ax.text(x, y, z,  num, color='red', ha='center', va='center', ma='center')

            COLOR = None
            for xyz_v in XYZ:
                xyz, v = xyz_v
                if v == 1:
                    COLOR = 'b'
                else:
                    COLOR = 'r'

                plt.plot(*xyz, COLOR, linewidth = 2 * linewidth)

            if eEFD is not None:
                if COLOR is None or COLOR == 'r':
                    COLOR = 'b'
                elif COLOR == 'b':
                    COLOR = 'r'
                else:
                    raise Exception()
                for xyz in eEFD:
                    plt.plot(*xyz, COLOR, linewidth=2 * linewidth)

            #------ We plot the involved trace-element-frame ----------------------------

            for ___ in Trace_Element_Frames:
                for lines in ___:
                    plt.plot(*lines, 'k', linewidth = 1.25*linewidth)

            for ___ in Trace_dof_xyz:
                for lines in ___:
                    plt.plot(*lines, '.', color='darkmagenta', linewidth = 3*linewidth)

            # --------- title ------------------------------------------------------------------
            if title is None:
                plt.title(f"Connection-(AD)1Trace-dof#{i}.")
            else:
                plt.title(title)

            if saveto is None:
                plt.show()
            else:
                plt.savefig(saveto, bbox_inches='tight')

            plt.close()
        # ======================================================================================

    def connection_through_around_edge_dof(self,
        i, T, C, tf, ef,
        density=5, zoom=0.9, linewidth=0.75, checking_mode=False, title=None,
        saveto=None,
        ):
        """"""
        #----------- we first do some checks -------------------------------------------------------
        assert self._sf_.whether.hybrid, f"the form is not hybrid, this method makes no sense."

        assert T.__class__.__name__ == 'EWC_SparseMatrix', f"T must be a trace matrix."
        assert C.__class__.__name__ == 'EWC_SparseMatrix', f"C must be a complement matrix."

        mesh = self._sf_.mesh

        if hasattr(tf, 'prime'): tf = tf.prime
        assert tf.__class__.__name__ == '_3dCSCG_1Trace'
        assert ef.__class__.__name__ == '_3dCSCG_1Edge'

        assert self._sf_.mesh == ef.mesh
        assert self._sf_.mesh == tf.mesh
        assert self._sf_.space == ef.space
        assert self._sf_.space == tf.space

        sGM = self._sf_.numbering.gathering
        tGM = tf.numbering.gathering
        eGM = ef.numbering.gathering
        ST = tf.matrices.selective.T

        # ------- the plot data of this edge dof -----------------------------------------------
        eDOF_PD = ef.special.generate_plot_data_for_dof(i, density=density, zoom=zoom) # in all cores

        #---------- we find the edge-element this dof is on ------------------------------------
        edge_element, local_index = ef.dofs.do.find.edge_element_and_local_index_of_dof(i)

        #-------- find the 1-sf, 1-tf dofs around this edge dof in all cores ----------------------------
        sf_dofs_numbers = list()
        tf_dofs_numbers = list()
        if edge_element in self._sf_.mesh.edge.elements:
            EE = self._sf_.mesh.edge.elements[edge_element]
            positions = EE.positions
            for pos in positions:
                if pos[0].isnumeric():
                    mesh_element = int(pos[:-2])
                    if mesh_element in mesh.elements:
                        corner_edge = pos[-2:]
                        sf_dof_local_index = self._sf_.numbering.do.find.\
                            local_dofs_on_element_corner_edge(corner_edge)[local_index]
                        sf_dof_number = sGM[mesh_element][sf_dof_local_index]
                        sf_dofs_numbers.append([mesh_element, sf_dof_number, sf_dof_local_index])

                        STm = ST[mesh_element]
                        STm_row = STm[sf_dof_local_index]
                        STm_row = STm_row.tocsr()
                        indices = STm_row.indices
                        assert len(indices) == 2, f"must be the case."
                        tf_dof_numbers = tGM[mesh_element][indices]
                        for tf_dof_local_index, tf_dof_number in zip(indices, tf_dof_numbers):
                            tf_dofs_numbers.append([tf_dof_number, tf_dof_local_index])
                else:
                    pass
        ___ = COMM.allgather(sf_dofs_numbers)
        sf_dofs_numbers = list()
        for _ in ___:
            sf_dofs_numbers.extend(_)

        ___ = COMM.allgather(tf_dofs_numbers)
        tf_dofs_numbers = list()
        for __ in ___:
            for _ in __:
                if _ not in tf_dofs_numbers:
                    tf_dofs_numbers.append(_)

        #-------- make the plot data of the 1-sf dofs only in master core --------------------------
        sDOF_PD = list()
        for _ in sf_dofs_numbers:
            mesh_element, sf_dof_number, sf_dof_local_index = _
            dof = self._dofs_[sf_dof_number]
            sDOF_pd = dof.do.generate_plot_data_from_element(
                mesh_element, density=2*density, zoom=zoom)

            efd = sDOF_pd['EFD']
            xyz = sDOF_pd['DOF_COO']

            if efd is not None:
                assert xyz is not None
                sDOF_PD.append([mesh_element,
                                sf_dof_number,
                                sf_dof_local_index,
                                efd,
                                xyz])

        sDOF_PD = COMM.gather(sDOF_PD, root=MASTER_RANK)
        if RANK == MASTER_RANK:
            ___ = list()
            for _ in sDOF_PD:
                ___.extend(_)
            sDOF_PD = ___

        #-------- make the plot data of the 1-tf dofs only in master core ------------------
        if RANK == MASTER_RANK:
            tDOF_PD = list()

        for _ in tf_dofs_numbers:
            tf_dof_number, tf_dof_local_index = _
            TEi_PD = tf.dofs[tf_dof_number].do.generate_plot_data(density=2*density, zoom=zoom)

            TEi_PD = COMM.gather(TEi_PD, root=MASTER_RANK)
            if RANK == MASTER_RANK:
                ___ = list()
                for _ in TEi_PD:
                    if _ is not None:
                        ___.append(_)
                TEi_PD = ___

                Trace_Element_Frames = list()
                Trace_dof_xyz = list()
                for _ in TEi_PD:
                    # noinspection PyTypeChecker
                    Trace_Element_Frames.append(_['Trace-Element-Frame'])
                    # noinspection PyTypeChecker
                    Trace_dof_xyz.append(_['DOF-xyz'])

                tDOF_PD.append([tf_dof_number,
                                tf_dof_local_index,
                                Trace_Element_Frames,
                                Trace_dof_xyz])

        #--------- do the checking ------------------------------------------------------
        if RANK == MASTER_RANK:
            involved_tf_dof = list()
            involved_sf_dofs = list()
            for tDOF_pd in tDOF_PD:
                involved_tf_dof.append(tDOF_pd[0]) # involved mesh-elements
            for sDOF_pd in sDOF_PD:
                involved_sf_dofs.append(sDOF_pd[1]) # involved mesh-elements
        else:
            involved_tf_dof = None
        involved_tf_dof = COMM.bcast(involved_tf_dof, root=MASTER_RANK)

        CONNECT = dict()
        involved_tf_dof = set(involved_tf_dof)
        for i_td in involved_tf_dof:
            t_elements_indices = tGM.do.find.elements_and_local_indices_of_dof(i_td)
            if t_elements_indices is not None:
                elements, indices = t_elements_indices
                for element, ind in zip(elements, indices):
                    T_row = T[element][ind].tocsr()
                    C_row = C[element][ind].tocsr()

                    if T_row.nnz > 0: # we find connection of sf dof for this trace dof.
                        indices = T_row.indices
                        assert len(indices) == 1, f"must connect only one sf dof."
                        if i_td not in CONNECT:
                            CONNECT[i_td] = list()

                        CONNECT[i_td].append(('sf', sGM[element][indices[0]], T_row[0,indices[0]]))
                    else:
                        pass

                    if C_row.nnz > 0: # we find connection of sf dof for this trace dof.
                        indices = C_row.indices
                        assert len(indices) == 1, f"must connect only one ef dof."
                        if i_td not in CONNECT:
                            CONNECT[i_td] = list()
                        assert eGM[element][indices[0]] == i, f"if find a ef dof, must be `i`."
                        CONNECT[i_td].append(('ef', i, C_row[0,indices[0]]))

        CONNECT = COMM.gather(CONNECT, root=MASTER_RANK)

        if RANK == MASTER_RANK:
            ___ = dict()
            for CON in CONNECT:
                for t_dof in CON:
                    if t_dof not in ___:
                        ___[t_dof] = list()
                    for _ in CON[t_dof]:
                        ___[t_dof].append(_)

            CONNECT = ___ # keys are trace dof number, values are the dofs the trace dofs connect.

            CT = Counter()

            # now we re-name these dofs
            DOF_dict = dict()
            for t_dof in CONNECT:
                for c_dof in CONNECT[t_dof]:
                    indicator = c_dof[0] + str(int(c_dof[1]))
                    if c_dof[0] == 'sf':
                        assert c_dof[1] in involved_sf_dofs, f'must be!'
                    elif c_dof[0] == 'ef':
                        assert c_dof[1] == i
                    if indicator not in DOF_dict:
                        DOF_dict[indicator] = next(CT)

            A = np.zeros((len(CONNECT), len(DOF_dict)))

            sf_dof_connecting_to_edge_dof = None
            for I, t_dof in enumerate(CONNECT):
                LEN = len(CONNECT[t_dof])
                for c_dof in CONNECT[t_dof]:
                    indicator = c_dof[0] + str(int(c_dof[1]))
                    value = c_dof[2]
                    j = DOF_dict[indicator]
                    A[I, j] = value

                if LEN == 2:
                    assert np.sum(A[I,:]) == 0, f"we must have found a +1 and a -1."

                    id0 = CONNECT[t_dof][0][0]
                    id1 = CONNECT[t_dof][1][0]

                    if id0 == 'sf' and id1 == 'ef':
                        sf_dof_connecting_to_edge_dof = CONNECT[t_dof][0][1]
                    elif id0 == 'ef' and id1 == 'sf':
                        sf_dof_connecting_to_edge_dof = CONNECT[t_dof][1][1]
                    else:
                        pass
                else:
                    assert LEN == 1

            if A.shape[0] > 0: # some edge dof may not be a part of the connection.
                A = rsmat(A)
                A = A[-1,:]
                if np.all(A == 0):
                    warnings.warn(f"Singularity find around edge dof#{i}, "
                                  f"do not forget applying BC later on",
                                  HybridSingularityWarning)

        #------------- make the plot ----------------------------------------------------------
        if RANK == MASTER_RANK and checking_mode is False:

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

            # ------ plot 1-sf dofs -----------------------------------------------------------
            SF_DOF_COLOR = ''

            for sDOF_pd in sDOF_PD:

                efd = sDOF_pd[3]
                if 'xLines_x' in efd:
                    X, Y, Z = efd['xLines_x'], efd['xLines_y'], efd['xLines_z']
                    for x, y, z in zip(X, Y, Z): plt.plot(x, y, z, 'lightgray', linewidth=linewidth)
                if 'yLines_x' in efd:
                    X, Y, Z = efd['yLines_x'], efd['yLines_y'], efd['yLines_z']
                    for x, y, z in zip(X, Y, Z): plt.plot(x, y, z, 'lightgray', linewidth=linewidth)
                if 'zLines_x' in efd:
                    X, Y, Z = efd['zLines_x'], efd['zLines_y'], efd['zLines_z']
                    for x, y, z in zip(X, Y, Z): plt.plot(x, y, z, 'lightgray', linewidth=linewidth)

                X, Y, Z = efd['xLines_x_B'], efd['xLines_y_B'], efd['xLines_z_B']
                for x, y, z in zip(X, Y, Z): plt.plot(x, y, z, 'gray', linewidth=linewidth)
                X, Y, Z = efd['yLines_x_B'], efd['yLines_y_B'], efd['yLines_z_B']
                for x, y, z in zip(X, Y, Z): plt.plot(x, y, z, 'gray', linewidth=linewidth)
                X, Y, Z = efd['zLines_x_B'], efd['zLines_y_B'], efd['zLines_z_B']
                for x, y, z in zip(X, Y, Z): plt.plot(x, y, z, 'gray', linewidth=linewidth)

                x, y, z = efd['center']['coordinate']
                x, y, z = x[0], y[0], z[0]
                num = efd['center']['number']
                ax.text(x, y, z,  num, color='red', ha='center', va='center', ma='center')

                xyz = sDOF_pd[4]

                if sf_dof_connecting_to_edge_dof is not None:

                    sf_dof_number = sDOF_pd[1]
                    if sf_dof_number == sf_dof_connecting_to_edge_dof:
                        plt.plot(*xyz, color='k', linewidth=2 * linewidth)
                        SF_DOF_COLOR = 'k'
                    else:
                        plt.plot(*xyz, linewidth=2 * linewidth)
                else:
                    plt.plot(*xyz, linewidth=2 * linewidth)

            # ------ plot 1-tf dofs -------------------------------------------------------------
            for tDOF_pd in tDOF_PD:
                Trace_Element_Frames = tDOF_pd[2]
                for ___ in Trace_Element_Frames:
                    for lines in ___:
                        plt.plot(*lines, 'gray', linewidth = 1.25*linewidth)

                Trace_dof_xyz = tDOF_pd[3]
                Trace_dof_number = tDOF_pd[0]
                for ___ in Trace_dof_xyz:
                    for lines in ___:
                        if Trace_dof_number in CONNECT:
                            plt.plot(*lines, '--', color='b', linewidth = 1.5*linewidth)
                        else:
                            plt.plot(*lines, '--', color='r', linewidth=1.5 * linewidth)

            # -- plot 1-edge form dof ----------------------------------------------------------
            for xyz in eDOF_PD:
                if SF_DOF_COLOR == 'k':
                    plt.plot(*xyz, color='k', linewidth=2 * linewidth)
                else: # this edge dof is not used, so skip plotting it.
                    pass

            # --------- title ------------------------------------------------------------------
            if title is None:
                plt.title(f"Connection-around-1-edge-dof#{i}.")
            else:
                plt.title(title)

            if saveto is None:
                plt.show()
            else:
                plt.savefig(saveto, bbox_inches='tight')

            plt.close()





if __name__ == '__main__':
    # mpiexec -n 5 python objects\CSCG\_3d\forms\standard\base\dofs\visualize\matplot\_1sf.py
    from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller

    elements = [2,2,2]
    # mesh = MeshGenerator('crazy_periodic', c=0.1)(elements)
    mesh = MeshGenerator('crazy', c=0.1)(elements)
    # Dirichlet_boundaries = []
    # Dirichlet_boundaries = ['North', 'South', 'West', 'East', 'Back', 'Front',] #
    # mesh = MeshGenerator('bridge_arch_cracked')(elements)
    space = SpaceInvoker('polynomials')([2, 2, 2])
    FC = FormCaller(mesh, space)

    f1 = FC('1-f', is_hybrid=True)
    t1 = FC('1-adt')
    e1 = FC('1-e')

    T = t1.matrices.trace
    C = e1.matrices.complement

    T, C = f1.special.___PRIVATE_overcoming_hybrid_singularity___(T, C)[:2]

    # f1.dofs.visualize.matplot.connection_through_around_edge_dof(55, T, C, t1, e1)


    # for i in range(t1.prime.numbering.gathering.GLOBAL_num_dofs):
    #     f1.dofs.visualize.matplot.connection_through_around_edge_dof(
    #         i, T, C, t1, e1)

    #
    for i in range(t1.prime.numbering.gathering.global_num_dofs):
        f1.dofs.visualize.matplot.connection_through_trace_dof(i, T, C, t1, e1, checking_mode=False)