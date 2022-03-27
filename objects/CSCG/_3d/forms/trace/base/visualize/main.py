


import tecplot as tp
from screws.freeze.main import FrozenOnly
from root.config.main import *
from tecplot.constant import PlotType
from objects.CSCG._3d.forms.trace.base.visualize.matplot import _3dCSCG_trace_form_Matplot



class _3dCSCG_Trace_Visualize(FrozenOnly):
    """The visualization property/component of standard forms."""
    def __init__(self, tf):
        self._tf_ = tf
        self._defaultPlot_ = 'matplot'
        self._matplot_ = _3dCSCG_trace_form_Matplot(tf)
        self._freeze_self_()

    def __call__(self, **kwargs):
        """When this object is called, we call the default visualizing method: ``tecplot``."""
        return getattr(self, self._defaultPlot_)(**kwargs)


    @property
    def matplot(self):
        """will visualize the trace form using ``matplotlib``."""
        return self._matplot_



    def tecplot(self, **kwargs):
        """Visualize the trace form using ``tecplot``."""
        getattr(self, f"_tecplot_{self._tf_.k}Trace_")(**kwargs)


    def _tecplot_2Trace_(self, numOfSamples=40000, port=7600):
        """Tecplot for 2-trace-form."""
        mesh = self._tf_.mesh
        numOfTotalTraceElements = mesh.trace.elements.GLOBAL_num
        density = int(np.ceil((numOfSamples / numOfTotalTraceElements) ** (1/2)))
        rst = [np.linspace(-1, 1, density) for _ in range(self._tf_.ndim)]
        xyz, v = self._tf_.reconstruct(*rst)
        # Now, we gather xyz & v from all cores into Master Core, store in XYZ & V ...
        xyz = cOmm.gather(xyz, root=mAster_rank)
        v = cOmm.gather(v, root=mAster_rank)
        if rAnk == mAster_rank:
            _ = dict()
            for xyz_i in xyz:
                _.update(xyz_i)
            xyz = _
            _ = dict()
            for v_i in v:
                _.update(v_i)
            v = _
        # Now, do the plot in the master core...
        if rAnk == mAster_rank:
            tp.session.connect(port=port)
            tp.new_layout()
            page = tp.active_page()
            frame = page.active_frame()
            page.name = '2traceFm: ' + self._tf_.standard_properties.name
            frame.name = '2traceFm: ' + self._tf_.standard_properties.name
            dataset = frame.create_dataset('2Trace')
            dataset.add_variable('x')
            dataset.add_variable('y')
            dataset.add_variable('z')
            dataset.add_variable('v')
            for key in v:
                # this makes the process slow. But anyway, it works. For now, we make keep it like this.
                size = np.shape(v[key][0])
                zone = dataset.add_ordered_zone('TraceElement_'+str(key), size)
                for j in range(3):
                    zone.values('xyz'[j])[:] = xyz[key][j].ravel('F')
                zone.values('v')[:] = v[key][0].ravel('F')
            frame.plot_type = getattr(PlotType, 'Cartesian3D')
            plot = frame.plot()
            plot.show_shade = False
            plot.show_edge = True
            plot.show_mesh = False
            plot.show_contour = True
            plot.use_translucency  = False
            plot.use_lighting_effect = False
            cont = plot.contour(0)
            cont.colormap_name = 'Diverging - Blue/Red'
            frame.plot().activate()