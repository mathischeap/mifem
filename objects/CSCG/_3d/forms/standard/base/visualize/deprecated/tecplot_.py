# -*- coding: utf-8 -*-



from root.config.main import *
import tecplot as tp
from tecplot.constant import PlotType
from screws.freeze.main import FrozenOnly



class _3dCSCG_standard_form_Tecplot(FrozenOnly):
    """"""

    def __init__(self, sf):
        """ """
        assert '3dCSCG_standard_form' in sf.standard_properties.tags
        self._sf_ = sf

        self._freeze_self_()


    def __call__(self, **kwargs):
        """Visualize the standard form using ``tecplot``."""
        getattr(self, f"_tecplot_{self._sf_.k}Form_")(**kwargs)


    def _tecplot_3Form_(self, numOfSamples=60000, port=7600):
        """ """
        mesh = self._sf_.mesh
        density = int(np.ceil((numOfSamples / mesh.elements.GLOBAL_num) ** (1/3)))
        rst = [np.linspace(-1, 1, density) for _ in range(self._sf_.ndim)]
        xyz, v = self._sf_.reconstruct(*rst)
        # Now, we gather xyz & v from all cores into Master Core, store in XYZ & V ...
        X = Y = Z = V = 0
        if rAnk == mAster_rank:
            X = [None for _ in range(mesh.elements.GLOBAL_num)]
            Y = [None for _ in range(mesh.elements.GLOBAL_num)]
            Z = [None for _ in range(mesh.elements.GLOBAL_num)]
            V = [None for _ in range(mesh.elements.GLOBAL_num)]
            for j in mesh.elements.indices:
                X[j] = xyz[j][0]
                Y[j] = xyz[j][1]
                Z[j] = xyz[j][2]
                V[j] = v[j][0]
            for i in sLave_ranks:
                xyz, v = cOmm.recv(source=i, tag=113)
                for j in xyz:
                    X[j] = xyz[j][0]
                    Y[j] = xyz[j][1]
                    Z[j] = xyz[j][2]
                    V[j] = v[j][0]
            del xyz, v
        else:
            cOmm.send([xyz, v], dest=mAster_rank, tag=113)
            del xyz, v

        # Now, we reshape the XYZ and V for tecplot and do the plot ...
        if rAnk == mAster_rank:
            X, Y, Z, V = mesh.do.regionwsie_stack(X, Y, Z, V)
            tp.session.connect(port=port)
            tp.new_layout()
            page = tp.active_page()
            frame = page.active_frame()
            page.name = '3stddFm: ' + self._sf_.standard_properties.name
            frame.name = '3stddFm: ' + self._sf_.standard_properties.name
            dataset = frame.create_dataset('3stddFm')
            for i in range(3): dataset.add_variable('xyz'[i])
            dataset.add_variable('v')
            for Rn in mesh.domain.regions.names:
                size = np.shape(V[Rn])
                zone = dataset.add_ordered_zone('Region_' + Rn, size)
                zone.values('x')[:] = X[Rn].ravel('F')
                zone.values('y')[:] = Y[Rn].ravel('F')
                zone.values('z')[:] = Z[Rn].ravel('F')
                zone.values('v')[:] = V[Rn].ravel('F')
            frame.plot_type = getattr(PlotType, 'Cartesian3D')
            plot = frame.plot()
            plot.show_shade = False
            plot.show_edge = True
            plot.show_mesh = False
            for j in range(mesh._num_total_elements_):
                surfaces = plot.fieldmap(j).surfaces
                surfaces.surfaces_to_plot = True
            plot.use_translucency  = False
            plot.use_lighting_effect = False
            plot.show_contour = True
            cont = plot.contour(0)
            cont.colormap_name = 'Diverging - Blue/Red'
            plot.view.fit()
            frame.plot().activate()



    def _tecplot_2Form_(self, numOfSamples=50000, port=7600):
        mesh = self._sf_.mesh
        density = int(np.ceil((numOfSamples / mesh.elements.GLOBAL_num) ** (1/3)))
        rst = [np.linspace(-1, 1, density) for _ in range(self._sf_.ndim)]
        xyz, v = self._sf_.reconstruct(*rst)
        # Now, we gather xyz & v from all cores into Master Core, store in XYZ & V ...
        X = Y = Z = Vx= Vy = Vz = 0
        if rAnk == mAster_rank:
            X = [None for _ in range(mesh.elements.GLOBAL_num)]
            Y = [None for _ in range(mesh.elements.GLOBAL_num)]
            Z = [None for _ in range(mesh.elements.GLOBAL_num)]
            Vx = [None for _ in range(mesh.elements.GLOBAL_num)]
            Vy = [None for _ in range(mesh.elements.GLOBAL_num)]
            Vz = [None for _ in range(mesh.elements.GLOBAL_num)]
            for j in mesh.elements.indices:
                X[j] = xyz[j][0]
                Y[j] = xyz[j][1]
                Z[j] = xyz[j][2]
                Vx[j] = v[j][0]
                Vy[j] = v[j][1]
                Vz[j] = v[j][2]
            for i in sLave_ranks:
                xyz, v = cOmm.recv(source=i, tag=112)
                for j in xyz:
                    X[j] = xyz[j][0]
                    Y[j] = xyz[j][1]
                    Z[j] = xyz[j][2]
                    Vx[j] = v[j][0]
                    Vy[j] = v[j][1]
                    Vz[j] = v[j][2]
            del xyz, v
        else:
            cOmm.send([xyz, v], dest=mAster_rank, tag=112)
            del xyz, v

        # Now, we reshape the XYZ and V for tecplot and do the plot ...
        if rAnk == mAster_rank:
            X, Y, Z, Vx, Vy, Vz = mesh.do.regionwsie_stack(X, Y, Z, Vx, Vy, Vz)
            tp.session.connect(port=port)
            tp.new_layout()
            page = tp.active_page()
            frame = page.active_frame()
            page.name = '2stddFm: ' + self._sf_.standard_properties.name
            frame.name = '2stddFm: ' + self._sf_.standard_properties.name
            dataset = frame.create_dataset('2stddFm')
            for i in range(3): dataset.add_variable('xyz'[i])
            dataset.add_variable('u')
            dataset.add_variable('v')
            dataset.add_variable('w')
            for Rn in mesh.domain.regions.names:
                size = np.shape(Vx[Rn])
                zone = dataset.add_ordered_zone('Region_' + Rn, size)
                zone.values('x')[:] = X[Rn].ravel('F')
                zone.values('y')[:] = Y[Rn].ravel('F')
                zone.values('z')[:] = Z[Rn].ravel('F')
                zone.values('u')[:] = Vx[Rn].ravel('F')
                zone.values('v')[:] = Vy[Rn].ravel('F')
                zone.values('w')[:] = Vz[Rn].ravel('F')
            frame.plot_type = getattr(PlotType, 'Cartesian3D')
            plot = frame.plot()
            plot.show_shade = False
            plot.show_edge = True
            plot.show_mesh = False
            for j in range(mesh._num_total_elements_):
                surfaces = plot.fieldmap(j).surfaces
                surfaces.surfaces_to_plot = True
            plot.use_translucency  = False
            plot.use_lighting_effect = False
            plot.show_contour = True
            cont = plot.contour(0)
            cont.colormap_name = 'Diverging - Blue/Red'
            plot.view.fit()
            frame.plot().activate()



    def _tecplot_1Form_(self, numOfSamples=50000, port=7600):
        mesh = self._sf_.mesh
        density = int(np.ceil((numOfSamples / mesh.elements.GLOBAL_num) ** (1/3)))
        rst = [np.linspace(-1, 1, density) for _ in range(self._sf_.ndim)]
        xyz, v = self._sf_.reconstruct(*rst)
        # Now, we gather xyz & v from all cores into Master Core, store in XYZ & V ...
        X = Y = Z = Vx= Vy = Vz = 0
        if rAnk == mAster_rank:
            X = [None for _ in range(mesh.elements.GLOBAL_num)]
            Y = [None for _ in range(mesh.elements.GLOBAL_num)]
            Z = [None for _ in range(mesh.elements.GLOBAL_num)]
            Vx = [None for _ in range(mesh.elements.GLOBAL_num)]
            Vy = [None for _ in range(mesh.elements.GLOBAL_num)]
            Vz = [None for _ in range(mesh.elements.GLOBAL_num)]
            for j in mesh.elements.indices:
                X[j] = xyz[j][0]
                Y[j] = xyz[j][1]
                Z[j] = xyz[j][2]
                Vx[j] = v[j][0]
                Vy[j] = v[j][1]
                Vz[j] = v[j][2]
            for i in sLave_ranks:
                xyz, v = cOmm.recv(source=i, tag=111)
                for j in xyz:
                    X[j] = xyz[j][0]
                    Y[j] = xyz[j][1]
                    Z[j] = xyz[j][2]
                    Vx[j] = v[j][0]
                    Vy[j] = v[j][1]
                    Vz[j] = v[j][2]
            del xyz, v
        else:
            cOmm.send([xyz, v], dest=mAster_rank, tag=111)
            del xyz, v

        # Now, we reshape the XYZ and V for tecplot and do the plot ...
        if rAnk == mAster_rank:
            X, Y, Z, Vx, Vy, Vz = mesh.do.regionwsie_stack(X, Y, Z, Vx, Vy, Vz)
            tp.session.connect(port=port)
            tp.new_layout()
            page = tp.active_page()
            frame = page.active_frame()
            page.name = '1stddFm: ' + self._sf_.standard_properties.name
            frame.name = '1stddFm: ' + self._sf_.standard_properties.name
            dataset = frame.create_dataset('1stddFm')
            for i in range(3): dataset.add_variable('xyz'[i])
            dataset.add_variable('u')
            dataset.add_variable('v')
            dataset.add_variable('w')
            for Rn in mesh.domain.regions.names:
                size = np.shape(Vx[Rn])
                zone = dataset.add_ordered_zone('Region_' + Rn, size)
                zone.values('x')[:] = X[Rn].ravel('F')
                zone.values('y')[:] = Y[Rn].ravel('F')
                zone.values('z')[:] = Z[Rn].ravel('F')
                zone.values('u')[:] = Vx[Rn].ravel('F')
                zone.values('v')[:] = Vy[Rn].ravel('F')
                zone.values('w')[:] = Vz[Rn].ravel('F')
            frame.plot_type = getattr(PlotType, 'Cartesian3D')
            plot = frame.plot()
            plot.show_shade = False
            plot.show_edge = True
            plot.show_mesh = False
            for j in range(mesh._num_total_elements_):
                surfaces = plot.fieldmap(j).surfaces
                surfaces.surfaces_to_plot = True
            plot.use_translucency  = False
            plot.use_lighting_effect = False
            plot.show_contour = True
            cont = plot.contour(0)
            cont.colormap_name = 'Diverging - Blue/Red'
            plot.view.fit()
            frame.plot().activate()



    def _tecplot_0Form_(self, numOfSamples=60000, port=7600):
        """ """
        mesh = self._sf_.mesh
        density = int(np.ceil((numOfSamples / mesh.elements.GLOBAL_num) ** (1/3)))
        rst = [np.linspace(-1, 1, density) for _ in range(self._sf_.ndim)]
        xyz, v = self._sf_.reconstruct(*rst)
        # Now, we gather xyz & v from all cores into Master Core, store in XYZ & V ...
        X = Y = Z = V = 0
        if rAnk == mAster_rank:
            X = [None for _ in range(mesh.elements.GLOBAL_num)]
            Y = [None for _ in range(mesh.elements.GLOBAL_num)]
            Z = [None for _ in range(mesh.elements.GLOBAL_num)]
            V = [None for _ in range(mesh.elements.GLOBAL_num)]
            for j in mesh.elements.indices:
                X[j] = xyz[j][0]
                Y[j] = xyz[j][1]
                Z[j] = xyz[j][2]
                V[j] = v[j][0]
            for i in sLave_ranks:
                xyz, v = cOmm.recv(source=i, tag=110)
                for j in xyz:
                    X[j] = xyz[j][0]
                    Y[j] = xyz[j][1]
                    Z[j] = xyz[j][2]
                    V[j] = v[j][0]
            del xyz, v
        else:
            cOmm.send([xyz, v], dest=mAster_rank, tag=110)
            del xyz, v

        # Now, we reshape the XYZ and V for tecplot and do the plot ...
        if rAnk == mAster_rank:
            X, Y, Z, V = mesh.do.regionwsie_stack(X, Y, Z, V)
            tp.session.connect(port=port)
            tp.new_layout()
            page = tp.active_page()
            frame = page.active_frame()
            page.name = '0stddFm: ' + self._sf_.standard_properties.name
            frame.name = '0stddFm: ' + self._sf_.standard_properties.name
            dataset = frame.create_dataset('0stddFm')
            for i in range(3): dataset.add_variable('xyz'[i])
            dataset.add_variable('v')
            for Rn in mesh.domain.regions.names:
                size = np.shape(V[Rn])
                zone = dataset.add_ordered_zone('Region_' + Rn, size)
                zone.values('x')[:] = X[Rn].ravel('F')
                zone.values('y')[:] = Y[Rn].ravel('F')
                zone.values('z')[:] = Z[Rn].ravel('F')
                zone.values('v')[:] = V[Rn].ravel('F')
            frame.plot_type = getattr(PlotType, 'Cartesian3D')
            plot = frame.plot()
            plot.show_shade = False
            plot.show_edge = True
            plot.show_mesh = False
            for j in range(mesh._num_total_elements_):
                surfaces = plot.fieldmap(j).surfaces
                surfaces.surfaces_to_plot = True
            plot.use_translucency  = False
            plot.use_lighting_effect = False
            plot.show_contour = True
            cont = plot.contour(0)
            cont.colormap_name = 'Diverging - Blue/Red'
            plot.view.fit()
            frame.plot().activate()



    def lambda_2_vortex(self, numOfSamples=60000, port=7600):
        """Plot the lambda_2 vortex detection.

        See [on the identification of a vortex] by Jeong and Hussain.
        """
        assert self._sf_.k in (1,2), f"lambda_2 vortex valid for 1- and 2-forms only."
        mesh = self._sf_.mesh
        density = int(np.ceil((numOfSamples / mesh.elements.GLOBAL_num) ** (1/3))) + 1
        rst = np.linspace(-1, 1, density)
        rst = 0.5 * (rst[1:] + rst[:-1])
        ___ = self._sf_.special.vortex_detection.Q_and_lambda2(rst, rst, rst)
        xyz, v = ___[0], ___[2] # xyz and lambda_2

        # Now, we gather xyz & v from all cores into Master Core, store in XYZ & V ...
        X = Y = Z = V = 0
        if rAnk == mAster_rank:
            X = [None for _ in range(mesh.elements.GLOBAL_num)]
            Y = [None for _ in range(mesh.elements.GLOBAL_num)]
            Z = [None for _ in range(mesh.elements.GLOBAL_num)]
            V = [None for _ in range(mesh.elements.GLOBAL_num)]
            for j in mesh.elements.indices:
                X[j] = xyz[j][0]
                Y[j] = xyz[j][1]
                Z[j] = xyz[j][2]
                V[j] = v[j]
            for i in sLave_ranks:
                xyz, v = cOmm.recv(source=i, tag=115)
                for j in xyz:
                    X[j] = xyz[j][0]
                    Y[j] = xyz[j][1]
                    Z[j] = xyz[j][2]
                    V[j] = v[j]
            del xyz, v
        else:
            cOmm.send([xyz, v], dest=mAster_rank, tag=115)
            del xyz, v

        # Now, we reshape the XYZ and V for tecplot and do the plot ...
        if rAnk == mAster_rank:
            X, Y, Z, V = mesh.do.regionwsie_stack(X, Y, Z, V)
            tp.session.connect(port=port)
            tp.new_layout()
            page = tp.active_page()
            frame = page.active_frame()
            page.name = 'lambda_2_vortex: ' + self._sf_.standard_properties.name
            frame.name = 'lambda_2_vortex: ' + self._sf_.standard_properties.name
            dataset = frame.create_dataset('lambda_2_vortex')
            for i in range(3): dataset.add_variable('xyz'[i])
            dataset.add_variable('v')
            for Rn in mesh.domain.regions.names:
                size = np.shape(V[Rn])
                zone = dataset.add_ordered_zone('Region_' + Rn, size)
                zone.values('x')[:] = X[Rn].ravel('F')
                zone.values('y')[:] = Y[Rn].ravel('F')
                zone.values('z')[:] = Z[Rn].ravel('F')
                zone.values('v')[:] = V[Rn].ravel('F')
            frame.plot_type = getattr(PlotType, 'Cartesian3D')
            plot = frame.plot()
            plot.show_shade = False
            plot.show_edge = True
            plot.show_mesh = False
            for j in range(mesh._num_total_elements_):
                surfaces = plot.fieldmap(j).surfaces
                surfaces.surfaces_to_plot = True
            plot.use_translucency  = False
            plot.use_lighting_effect = False
            plot.show_contour = True
            cont = plot.contour(0)
            cont.colormap_name = 'Diverging - Blue/Red'
            plot.view.fit()
            frame.plot().activate()



    def Q_vortex(self, numOfSamples=60000, port=7600):
        """Plot the Q vortex detection.

        See [on the identification of a vortex] by Jeong and Hussain.
        """
        assert self._sf_.k in (1,2), f"Q vortex valid for 1- and 2-forms only."
        mesh = self._sf_.mesh
        density = int(np.ceil((numOfSamples / mesh.elements.GLOBAL_num) ** (1 / 3))) + 1
        rst = np.linspace(-1, 1, density)
        rst = 0.5 * (rst[1:] + rst[:-1])
        ___ = self._sf_.special.vortex_detection.Q_and_lambda2(rst, rst, rst)
        xyz, v = ___[0], ___[1] # xyz and Q

        # Now, we gather xyz & v from all cores into Master Core, store in XYZ & V ...
        X = Y = Z = V = 0
        if rAnk == mAster_rank:
            X = [None for _ in range(mesh.elements.GLOBAL_num)]
            Y = [None for _ in range(mesh.elements.GLOBAL_num)]
            Z = [None for _ in range(mesh.elements.GLOBAL_num)]
            V = [None for _ in range(mesh.elements.GLOBAL_num)]
            for j in mesh.elements.indices:
                X[j] = xyz[j][0]
                Y[j] = xyz[j][1]
                Z[j] = xyz[j][2]
                V[j] = v[j]
            for i in sLave_ranks:
                xyz, v = cOmm.recv(source=i, tag=112)
                for j in xyz:
                    X[j] = xyz[j][0]
                    Y[j] = xyz[j][1]
                    Z[j] = xyz[j][2]
                    V[j] = v[j]
            del xyz, v
        else:
            cOmm.send([xyz, v], dest=mAster_rank, tag=112)
            del xyz, v

        # Now, we reshape the XYZ and V for tecplot and do the plot ...
        if rAnk == mAster_rank:
            X, Y, Z, V = mesh.do.regionwsie_stack(X, Y, Z, V)
            tp.session.connect(port=port)
            tp.new_layout()
            page = tp.active_page()
            frame = page.active_frame()
            page.name = 'Q_vortex: ' + self._sf_.standard_properties.name
            frame.name = 'Q_vortex: ' + self._sf_.standard_properties.name
            dataset = frame.create_dataset('Q_vortex')
            for i in range(3): dataset.add_variable('xyz'[i])
            dataset.add_variable('v')
            for Rn in mesh.domain.regions.names:
                size = np.shape(V[Rn])
                zone = dataset.add_ordered_zone('Region_' + Rn, size)
                zone.values('x')[:] = X[Rn].ravel('F')
                zone.values('y')[:] = Y[Rn].ravel('F')
                zone.values('z')[:] = Z[Rn].ravel('F')
                zone.values('v')[:] = V[Rn].ravel('F')
            frame.plot_type = getattr(PlotType, 'Cartesian3D')
            plot = frame.plot()
            plot.show_shade = False
            plot.show_edge = True
            plot.show_mesh = False
            for j in range(mesh._num_total_elements_):
                surfaces = plot.fieldmap(j).surfaces
                surfaces.surfaces_to_plot = True
            plot.use_translucency  = False
            plot.use_lighting_effect = False
            plot.show_contour = True
            cont = plot.contour(0)
            cont.colormap_name = 'Diverging - Blue/Red'
            plot.view.fit()
            frame.plot().activate()