


import sys
if './' not in sys.path: sys.path.append('/')

from screws.frozen import FrozenOnly
from root.config import *

import matplotlib.pyplot as plt






class _3dCSCG_SF_DOF_VISUALIZE(FrozenOnly):
    """"""
    def __init__(self, dof):
        """"""
        self._dof_ = dof
        self._mesh_ = dof._sf_.mesh
        self._matplot_ = _3dCSCG_SF_DOF_VISUALIZE_matplot(dof)
        self._freeze_self_()

    def __call__(self, *args, **kwargs):
        return self.matplot(*args, **kwargs)

    @property
    def matplot(self):
        return self._matplot_





class _3dCSCG_SF_DOF_VISUALIZE_matplot(FrozenOnly):
    """"""
    def __init__(self, dof):
        """"""
        self._dof_ = dof
        self._mesh_ = dof._sf_.mesh
        self._sf_ = dof._sf_
        self._freeze_self_()

    def __call__(self, *args, **kwargs):
        """Default visualizer"""
        return getattr(self, f"___PRIVATE_matplot_dof_{self._sf_.k}form___")(*args, **kwargs)



    def ___PRIVATE_matplot_dof_0form___(self, *args, **kwargs):
        """We plot this dof of a standard 0-form."""
        if self._sf_.IS_hybrid:
            return self.___PRIVATE_matplot_dof_k_form_IS_hybrid__(*args, **kwargs)
        else:
            return self.___PRIVATE_matplot_dof_0form_IS_NOT_hybrid__(*args, **kwargs)

    def ___PRIVATE_matplot_dof_0form_IS_NOT_hybrid__(self, density=20, saveto=None, linewidth=0.6, title=None):
        """"""
        positions = self._dof_.positions
        EF = dict()
        for E_I in positions:
            E, I = E_I
            EF[E] = self._sf_.___PRIVATE_element_grid_data_generator_1___(E, density=density)
        GPs = self._dof_.GLOBAL_positions
        Element_Frames = cOmm.gather(EF, root=mAster_rank)
        if rAnk == mAster_rank:
            #------------ prepare figure -----------------------------------------------------------------
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
            ax.tick_params(labelsize=12)
            ax.set_xlabel(r'$x$', fontsize=15)
            ax.set_ylabel(r'$y$', fontsize=15)
            ax.set_zlabel(r'$z$', fontsize=15)
            #------ plot element frame -----------------------------------------------------------------
            Element_Frame = dict()
            for ef in Element_Frames:
                Element_Frame.update(ef)
            for e in Element_Frame:
                if 'xLines_x' in Element_Frame[e]:
                    X, Y, Z = Element_Frame[e]['xLines_x'], Element_Frame[e]['xLines_y'], Element_Frame[e]['xLines_z']
                    for x, y, z in zip(X, Y, Z):
                        plt.plot(x, y, z, 'gray', linewidth=linewidth)

                if 'yLines_x' in Element_Frame[e]:
                    X, Y, Z = Element_Frame[e]['yLines_x'], Element_Frame[e]['yLines_y'], Element_Frame[e]['yLines_z']
                    for x, y, z in zip(X, Y, Z):
                        plt.plot(x, y, z, 'gray', linewidth=linewidth)

                if 'zLines_x' in Element_Frame[e]:
                    X, Y, Z = Element_Frame[e]['zLines_x'], Element_Frame[e]['zLines_y'], Element_Frame[e]['zLines_z']
                    for x, y, z in zip(X, Y, Z):
                        plt.plot(x, y, z, 'gray', linewidth=linewidth)

                X, Y, Z = Element_Frame[e]['xLines_x_B'], Element_Frame[e]['xLines_y_B'], Element_Frame[e]['xLines_z_B']
                for x, y, z in zip(X, Y, Z):
                    plt.plot(x, y, z, 'green', linewidth=linewidth)
                X, Y, Z = Element_Frame[e]['yLines_x_B'], Element_Frame[e]['yLines_y_B'], Element_Frame[e]['yLines_z_B']
                for x, y, z in zip(X, Y, Z):
                    plt.plot(x, y, z, 'green', linewidth=linewidth)
                X, Y, Z = Element_Frame[e]['zLines_x_B'], Element_Frame[e]['zLines_y_B'], Element_Frame[e]['zLines_z_B']
                for x, y, z in zip(X, Y, Z):
                    plt.plot(x, y, z, 'green', linewidth=linewidth)

                x, y, z = Element_Frame[e]['center']['coordinate']
                x, y, z = x[0], y[0], z[0]
                num = Element_Frame[e]['center']['number']
                ax.text(x, y, z,  num, color='red', ha='center', va='center', ma='center')

            #------ plot the dof of 0-form ----------------------------------------------------------------
            pos = GPs[0]
            element, index = pos
            i, j, k = np.where(self._sf_.numbering.local[0]==index)
            i, j, k = i[0], j[0], k[0]
            nodes = self._sf_.space.nodes
            x, y, z = nodes[0][i], nodes[1][j], nodes[2][k]
            x = np.array([x])
            y = np.array([y])
            z = np.array([z])
            xyz = x, y, z
        else:
            element = None
            xyz = None
        xyz, element = cOmm.bcast([xyz, element], root=mAster_rank)

        if element in self._mesh_.elements:
            xyz = self._mesh_.elements[element].coordinate_transformation.mapping(*xyz)
        else:
            xyz = None

        xyz = cOmm.gather(xyz, root=mAster_rank)

        if rAnk == mAster_rank:
            for ___ in xyz:
                if ___ is not None:
                    x, y, z = ___
                    break

            ax.scatter(x[0], y[0], z[0], marker='s', color='b')

            #--------- title ------------------------------------------------------------------------------
            if title is None:
                plt.title(f"dof#{self._dof_._i_} of {self._sf_.k}-form: {self._sf_.standard_properties.name}.")
            else:
                plt.title(title)

            #---------- SAVE TO ---------------------------------------------------------------------------
            plt.tight_layout()
            if saveto is not None and saveto != '':
                plt.savefig(saveto, bbox_inches='tight')
                plt.close()
            else:
                plt.show()
            #================================================================================

            return fig

    def ___PRIVATE_matplot_dof_k_form_IS_hybrid__(self, *args, **kwargs):
        """"""
        f_kwargs = self._sf_.___define_parameters___['kwargs']
        f_kwargs['is_hybrid'] = False
        non_hybrid_form = self._sf_.__class__(self._sf_.mesh, self._sf_.space, **f_kwargs)
        dofs = non_hybrid_form.dofs
        hy_i = self._dof_._i_

        GPs = self._dof_.GLOBAL_positions
        assert len(GPs) == 1, f"trivial check!"

        pos = GPs[0]
        if pos in self._dof_.positions:
            Ele, index = pos
            nhy_i = non_hybrid_form.numbering.gathering[Ele][index]
        else:
            nhy_i = None

        nhy_i = cOmm.gather(nhy_i, root=mAster_rank)
        if rAnk == mAster_rank:
            I = 0
            for _ in nhy_i:
                if _ is not None:
                    I += 1
                    NHY_I = _
            assert I == 1, "only find one position."

            nhy_i = NHY_I
        else:
            pass

        nhy_i = cOmm.bcast(nhy_i, root=mAster_rank)
        DI = dofs[nhy_i]

        if 'title' not in kwargs:
            kwargs['title'] = f"dof#{hy_i} of {self._sf_.k}-form: {self._sf_.standard_properties.name}."
        DI.visualize.matplot(*args, **kwargs)



    def ___PRIVATE_matplot_dof_1form___(self, *args, **kwargs):
        """We plot this dof of a standard 0-form."""
        if self._sf_.IS_hybrid:
            return self.___PRIVATE_matplot_dof_k_form_IS_hybrid__(*args, **kwargs)
        else:
            return self.___PRIVATE_matplot_dof_1form_IS_NOT_hybrid__(*args, **kwargs)

    def ___PRIVATE_matplot_dof_1form_IS_NOT_hybrid__(self, density=20, saveto=None, linewidth=0.6, title=None):
        """"""
        positions = self._dof_.positions
        EF = dict()
        for E_I in positions:
            E, I = E_I
            EF[E] = self._sf_.___PRIVATE_element_grid_data_generator_1___(E, density=density)
        GPs = self._dof_.GLOBAL_positions
        Element_Frames = cOmm.gather(EF, root=mAster_rank)
        if rAnk == mAster_rank:
            # ------------ prepare figure -----------------------------------------------------------------
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
            # ------ plot element frame -----------------------------------------------------------------
            Element_Frame = dict()
            for ef in Element_Frames:
                Element_Frame.update(ef)
            for e in Element_Frame:
                if 'xLines_x' in Element_Frame[e]:
                    X, Y, Z = Element_Frame[e]['xLines_x'], Element_Frame[e]['xLines_y'], Element_Frame[e]['xLines_z']
                    for x, y, z in zip(X, Y, Z):
                        plt.plot(x, y, z, 'gray', linewidth=linewidth)

                if 'yLines_x' in Element_Frame[e]:
                    X, Y, Z = Element_Frame[e]['yLines_x'], Element_Frame[e]['yLines_y'], Element_Frame[e]['yLines_z']
                    for x, y, z in zip(X, Y, Z):
                        plt.plot(x, y, z, 'gray', linewidth=linewidth)

                if 'zLines_x' in Element_Frame[e]:
                    X, Y, Z = Element_Frame[e]['zLines_x'], Element_Frame[e]['zLines_y'], Element_Frame[e]['zLines_z']
                    for x, y, z in zip(X, Y, Z):
                        plt.plot(x, y, z, 'gray', linewidth=linewidth)

                X, Y, Z = Element_Frame[e]['xLines_x_B'], Element_Frame[e]['xLines_y_B'], Element_Frame[e]['xLines_z_B']
                for x, y, z in zip(X, Y, Z):
                    plt.plot(x, y, z, 'green', linewidth=linewidth)
                X, Y, Z = Element_Frame[e]['yLines_x_B'], Element_Frame[e]['yLines_y_B'], Element_Frame[e]['yLines_z_B']
                for x, y, z in zip(X, Y, Z):
                    plt.plot(x, y, z, 'green', linewidth=linewidth)
                X, Y, Z = Element_Frame[e]['zLines_x_B'], Element_Frame[e]['zLines_y_B'], Element_Frame[e]['zLines_z_B']
                for x, y, z in zip(X, Y, Z):
                    plt.plot(x, y, z, 'green', linewidth=linewidth)

                x, y, z = Element_Frame[e]['center']['coordinate']
                x, y, z = x[0], y[0], z[0]
                num = Element_Frame[e]['center']['number']
                ax.text(x, y, z,  num, color='red', ha='center', va='center', ma='center')

            # ------ plot the dof of 1-form ----------------------------------------------------------------
            pos = GPs[0]
            element, index = pos
            i, j, k = np.where(self._sf_.numbering.local[0] == index)
            if len(i) == 0:
                i, j, k = np.where(self._sf_.numbering.local[1] == index)
                if len(i) == 0:
                    i, j, k = np.where(self._sf_.numbering.local[2] == index)
                    assert len(i) != 0, f"Must have found the index."
                    IND = 2
                else:
                    IND = 1
            else:
                IND = 0
            i, j, k = i[0], j[0], k[0]
            nodes = self._sf_.space.nodes
            if IND == 0:
                x, y, z = np.linspace(nodes[0][i], nodes[0][i+1], density), \
                          nodes[1][j] * np.ones(density), \
                          nodes[2][k] * np.ones(density)
            elif IND == 1:
                x, y, z = nodes[0][i] * np.ones(density), \
                          np.linspace(nodes[1][j], nodes[1][j+1], density), \
                          nodes[2][k] * np.ones(density)
            elif IND == 2:
                x, y, z = nodes[0][i] * np.ones(density), \
                          nodes[1][j] * np.ones(density), \
                          np.linspace(nodes[2][k], nodes[2][k+1], density)
            else:
                raise Exception()

            xyz = x, y, z
        else:
            element = None
            xyz = None
        xyz, element = cOmm.bcast([xyz, element], root=mAster_rank)

        if element in self._mesh_.elements:
            xyz = self._mesh_.elements[element].coordinate_transformation.mapping(*xyz)
        else:
            xyz = None

        xyz = cOmm.gather(xyz, root=mAster_rank)

        if rAnk == mAster_rank:
            for ___ in xyz:
                if ___ is not None:
                    x, y, z = ___
                    break
            plt.plot(x, y, z, 'blue', linewidth=3*linewidth)
            # ax.scatter(x[0], y[0], z[0], marker='s', color='b')

            # --------- title ------------------------------------------------------------------------------
            if title is None:
                plt.title(f"dof#{self._dof_._i_} of {self._sf_.k}-form: {self._sf_.standard_properties.name}.")
            else:
                plt.title(title)
            # ---------- SAVE TO ---------------------------------------------------------------------------
            plt.tight_layout()
            if saveto is not None and saveto != '':
                plt.savefig(saveto, bbox_inches='tight')
                plt.close()
            else:
                plt.show()
            # ================================================================================

            return fig



    def ___PRIVATE_matplot_dof_2form___(self, *args, **kwargs):
        """We plot this dof of a standard 0-form."""
        if self._sf_.IS_hybrid:
            return self.___PRIVATE_matplot_dof_k_form_IS_hybrid__(*args, **kwargs)
        else:
            return self.___PRIVATE_matplot_dof_2form_IS_NOT_hybrid__(*args, **kwargs)

    def ___PRIVATE_matplot_dof_2form_IS_NOT_hybrid__(self, density=20, saveto=None, linewidth=0.6, title=None):
        """"""
        positions = self._dof_.positions
        EF = dict()
        for E_I in positions:
            E, I = E_I
            EF[E] = self._sf_.___PRIVATE_element_grid_data_generator_1___(E, density=density)
        GPs = self._dof_.GLOBAL_positions
        Element_Frames = cOmm.gather(EF, root=mAster_rank)
        if rAnk == mAster_rank:
            # ------------ prepare figure -----------------------------------------------------------------
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
            # ------ plot element frame -----------------------------------------------------------------
            Element_Frame = dict()
            for ef in Element_Frames:
                Element_Frame.update(ef)
            for e in Element_Frame:
                if 'xLines_x' in Element_Frame[e]:
                    X, Y, Z = Element_Frame[e]['xLines_x'], Element_Frame[e]['xLines_y'], Element_Frame[e]['xLines_z']
                    for x, y, z in zip(X, Y, Z):
                        plt.plot(x, y, z, 'gray', linewidth=linewidth)

                if 'yLines_x' in Element_Frame[e]:
                    X, Y, Z = Element_Frame[e]['yLines_x'], Element_Frame[e]['yLines_y'], Element_Frame[e]['yLines_z']
                    for x, y, z in zip(X, Y, Z):
                        plt.plot(x, y, z, 'gray', linewidth=linewidth)

                if 'zLines_x' in Element_Frame[e]:
                    X, Y, Z = Element_Frame[e]['zLines_x'], Element_Frame[e]['zLines_y'], Element_Frame[e]['zLines_z']
                    for x, y, z in zip(X, Y, Z):
                        plt.plot(x, y, z, 'gray', linewidth=linewidth)

                X, Y, Z = Element_Frame[e]['xLines_x_B'], Element_Frame[e]['xLines_y_B'], Element_Frame[e]['xLines_z_B']
                for x, y, z in zip(X, Y, Z):
                    plt.plot(x, y, z, 'green', linewidth=linewidth)
                X, Y, Z = Element_Frame[e]['yLines_x_B'], Element_Frame[e]['yLines_y_B'], Element_Frame[e]['yLines_z_B']
                for x, y, z in zip(X, Y, Z):
                    plt.plot(x, y, z, 'green', linewidth=linewidth)
                X, Y, Z = Element_Frame[e]['zLines_x_B'], Element_Frame[e]['zLines_y_B'], Element_Frame[e]['zLines_z_B']
                for x, y, z in zip(X, Y, Z):
                    plt.plot(x, y, z, 'green', linewidth=linewidth)

                x, y, z = Element_Frame[e]['center']['coordinate']
                x, y, z = x[0], y[0], z[0]
                num = Element_Frame[e]['center']['number']
                ax.text(x, y, z,  num, color='red', ha='center', va='center', ma='center')

            # ------ plot the dof of 2-form ----------------------------------------------------------------
            pos = GPs[0]
            element, index = pos
            i, j, k = np.where(self._sf_.numbering.local[0] == index)
            if len(i) == 0:
                i, j, k = np.where(self._sf_.numbering.local[1] == index)
                if len(i) == 0:
                    i, j, k = np.where(self._sf_.numbering.local[2] == index)
                    assert len(i) != 0, f"Must have found the index."
                    IND = 2
                else:
                    IND = 1
            else:
                IND = 0
            i, j, k = i[0], j[0], k[0]
            nodes = self._sf_.space.nodes #
            if IND == 0:
                x, y, z = nodes[0][i] * np.ones((density, density)), \
                          np.linspace(nodes[1][j], nodes[1][j+1], density), \
                          np.linspace(nodes[2][k], nodes[2][k+1], density)
                y, z = np.meshgrid(y, z, indexing='ij')
            elif IND == 1:
                x, y, z = np.linspace(nodes[0][i], nodes[0][i+1], density), \
                          nodes[1][j] * np.ones((density, density)), \
                          np.linspace(nodes[2][k], nodes[2][k+1], density)
                x, z = np.meshgrid(x, z, indexing='ij')
            elif IND == 2:
                x, y, z = np.linspace(nodes[0][i], nodes[0][i+1], density), \
                          np.linspace(nodes[1][j], nodes[1][j+1], density), \
                          nodes[2][k] * np.ones((density, density))
                x, y = np.meshgrid(x, y, indexing='ij')
            else:
                raise Exception()
            xyz = x, y, z
        else:
            element = None
            xyz = None
        xyz, element = cOmm.bcast([xyz, element], root=mAster_rank)

        if element in self._mesh_.elements:
            xyz = self._mesh_.elements[element].coordinate_transformation.mapping(*xyz)
        else:
            xyz = None

        xyz = cOmm.gather(xyz, root=mAster_rank)

        if rAnk == mAster_rank:
            for ___ in xyz:
                if ___ is not None:
                    x, y, z = ___
                    break

            ax.plot_surface(x, y, z, color=(0,0,1,0.5))

            # --------- title ------------------------------------------------------------------------------
            if title is None:
                plt.title(f"dof#{self._dof_._i_} of {self._sf_.k}-form: {self._sf_.standard_properties.name}.")
            else:
                plt.title(title)
            # ---------- SAVE TO ---------------------------------------------------------------------------
            plt.tight_layout()
            if saveto is not None and saveto != '':
                plt.savefig(saveto, bbox_inches='tight')
                plt.close()
            else:
                plt.show()
            # ================================================================================

            return fig



    def ___PRIVATE_matplot_dof_3form___(self, density=15, saveto=None, linewidth=0.6, title=None):
        """No matter it is hybrid or not."""
        positions = self._dof_.positions
        EF = dict()
        for E_I in positions:
            E, I = E_I
            EF[E] = self._sf_.___PRIVATE_element_grid_data_generator_1___(E, density=density)
        GPs = self._dof_.GLOBAL_positions
        assert len(GPs) == 1, f"trivial check."
        Element_Frames = cOmm.gather(EF, root=mAster_rank)
        if rAnk == mAster_rank:
            # ------------ prepare figure -----------------------------------------------------------------
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
            # ------ plot element frame -----------------------------------------------------------------
            Element_Frame = dict()
            for ef in Element_Frames:
                Element_Frame.update(ef)
            assert len(Element_Frame) == 1, f"trivial check!"
            for e in Element_Frame:
                if 'xLines_x' in Element_Frame[e]:
                    X, Y, Z = Element_Frame[e]['xLines_x'], Element_Frame[e]['xLines_y'], Element_Frame[e]['xLines_z']
                    for x, y, z in zip(X, Y, Z):
                        plt.plot(x, y, z, 'gray', linewidth=linewidth)

                if 'yLines_x' in Element_Frame[e]:
                    X, Y, Z = Element_Frame[e]['yLines_x'], Element_Frame[e]['yLines_y'], Element_Frame[e]['yLines_z']
                    for x, y, z in zip(X, Y, Z):
                        plt.plot(x, y, z, 'gray', linewidth=linewidth)

                if 'zLines_x' in Element_Frame[e]:
                    X, Y, Z = Element_Frame[e]['zLines_x'], Element_Frame[e]['zLines_y'], Element_Frame[e]['zLines_z']
                    for x, y, z in zip(X, Y, Z):
                        plt.plot(x, y, z, 'gray', linewidth=linewidth)

                X, Y, Z = Element_Frame[e]['xLines_x_B'], Element_Frame[e]['xLines_y_B'], Element_Frame[e]['xLines_z_B']
                for x, y, z in zip(X, Y, Z):
                    plt.plot(x, y, z, 'green', linewidth=linewidth)
                X, Y, Z = Element_Frame[e]['yLines_x_B'], Element_Frame[e]['yLines_y_B'], Element_Frame[e]['yLines_z_B']
                for x, y, z in zip(X, Y, Z):
                    plt.plot(x, y, z, 'green', linewidth=linewidth)
                X, Y, Z = Element_Frame[e]['zLines_x_B'], Element_Frame[e]['zLines_y_B'], Element_Frame[e]['zLines_z_B']
                for x, y, z in zip(X, Y, Z):
                    plt.plot(x, y, z, 'green', linewidth=linewidth)

                x, y, z = Element_Frame[e]['center']['coordinate']
                x, y, z = x[0], y[0], z[0]
                num = Element_Frame[e]['center']['number']
                ax.text(x, y, z,  num, color='red', ha='center', va='center', ma='center')

            # ------ plot the dof of 2-form ----------------------------------------------------------------
            pos = GPs[0]
            element, index = pos
            i, j, k = np.where(self._sf_.numbering.local[0] == index)
            i, j, k = i[0], j[0], k[0]
            nodes = self._sf_.space.nodes  #
            x, y, z = np.linspace(nodes[0][i], nodes[0][i + 1], density), \
                      np.linspace(nodes[1][j], nodes[1][j + 1], density), \
                      np.linspace(nodes[2][k], nodes[2][k + 1], density)
            xyz = np.meshgrid(x, y, z, indexing='ij')
        else:
            element = None
            xyz = None
        xyz, element = cOmm.bcast([xyz, element], root=mAster_rank)

        if element in self._mesh_.elements:
            xyz = self._mesh_.elements[element].coordinate_transformation.mapping(*xyz)
        else:
            xyz = None

        xyz = cOmm.gather(xyz, root=mAster_rank)

        if rAnk == mAster_rank:
            for ___ in xyz:
                if ___ is not None:
                    X, Y, Z = ___
                    break

            # North: x-
            x, y, z = X[0, :, :], Y[0, :, :], Z[0, :, :],
            ax.plot_surface(x, y, z, color=(0, 0, 1, 0.3))
            # South: x+
            x, y, z = X[-1, :, :], Y[-1, :, :], Z[-1, :, :],
            ax.plot_surface(x, y, z, color=(0, 0, 1, 0.3))
            # West: x-
            x, y, z = X[:, 0, :], Y[:, 0, :], Z[:, 0, :],
            ax.plot_surface(x, y, z, color=(0, 0, 1, 0.3))
            # East: x+
            x, y, z = X[:, -1, :], Y[:, -1, :], Z[:, -1, :],
            ax.plot_surface(x, y, z, color=(0, 0, 1, 0.3))
            # Back: x-
            x, y, z = X[:, :, 0], Y[:, :, 0], Z[:, :, 0],
            ax.plot_surface(x, y, z, color=(0, 0, 1, 0.3))
            # Front: x+
            x, y, z = X[:, :, -1], Y[:, :, -1], Z[:, :, -1],
            ax.plot_surface(x, y, z, color=(0, 0, 1, 0.3))

            # --------- title ------------------------------------------------------------------------------
            if title is None:
                plt.title(
                    f"dof#{self._dof_._i_} of {self._sf_.k}-form: {self._sf_.standard_properties.name}.")
            else:
                plt.title(title)
            # ---------- SAVE TO ---------------------------------------------------------------------------
            plt.tight_layout()
            if saveto is not None and saveto != '':
                plt.savefig(saveto, bbox_inches='tight')
                plt.close()
            else:
                plt.show()
            # ================================================================================

            return fig









if __name__ == '__main__':
    # mpiexec -n 6 python _3dCSCG\form\standard\dofs\dof\main.py
    from _3dCSCG.main import MeshGenerator, SpaceInvoker, FormCaller#, ExactSolutionSelector

    mesh = MeshGenerator('crazy', c=0.3)([3,3,3])
    space = SpaceInvoker('polynomials')([('Lobatto',3), ('Lobatto',3), ('Lobatto',3)])
    FC = FormCaller(mesh, space)

    f = FC('3-f', is_hybrid=True)

    # SPL = f0.numbering.sharing_physical_locations
    #
    # print(SPL)

    dofs = f.dofs
    for i in dofs:
        DI = dofs[i]
        DI.visualize()