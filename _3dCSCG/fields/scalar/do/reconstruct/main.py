from screws.freeze.inheriting.frozen_only import FrozenOnly
from _3dCSCG.fields.scalar.do.reconstruct.mesh_element.standard import OnMeshElement_Standard
from _3dCSCG.fields.scalar.do.reconstruct.trace_element.trace_element_wise import OnTraceElement_TraceElementWise
from _3dCSCG.fields.scalar.do.reconstruct.trace_element.boundary_wise import OnTraceElement_BoundaryWise

class _3dCSCG_Scalar_Do_Reconstruct(FrozenOnly):
    def __init__(self, sf):
        self._sf_ = sf
        self._on_mesh_element___for_standard_ = OnMeshElement_Standard(sf)
        self._on_trace_element___for_trace_element_wise_ = OnTraceElement_TraceElementWise(sf)
        self._on_trace_element___for_boundary_wise_ = OnTraceElement_BoundaryWise(sf)
        self._freeze_self_()


    def __call__(self, xi, eta, sigma, time=None, ravel=False, i=None, where=None):
        """

        :param xi:
        :param eta:
        :param sigma:
        :param time:
        :param ravel:
        :param i:
        :param where:
        :return:
        """
        ftype = self._sf_.ftype

        #------- deal with time --------------------------------------------------
        if time is None:
            pass
        else:
            self._sf_.current_time = time

        # we deal with default `where` input ---------------------------------------------------------------
        if where is None:
            if ftype == "standard":
                where = "mesh-element"
            elif ftype in ("boundary-wise", "trace-element-wise"):
                where = "trace-element"
            else:
                raise NotImplementedError()
        else:
            pass

        #--------------------------------------------------------------------------
        if where == 'mesh-element':
            if ftype == 'standard':
                return self._on_mesh_element___for_standard_(xi, eta, sigma, ravel, i)
            else:
                raise NotImplementedError()
        elif where == 'trace-element':
            if ftype == 'trace-element-wise':
                return self._on_trace_element___for_trace_element_wise_(xi, eta, sigma, ravel, i)
            elif ftype == 'boundary-wise':
                return self._on_trace_element___for_boundary_wise_(xi, eta, sigma, ravel, i)
            else:
                raise NotImplementedError()
        else:
            raise NotImplementedError()
