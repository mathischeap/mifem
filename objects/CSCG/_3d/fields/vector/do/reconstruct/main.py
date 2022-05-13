# -*- coding: utf-8 -*-
from screws.freeze.base import FrozenOnly
from objects.CSCG._3d.fields.vector.do.reconstruct.mesh_element.standard import OnMeshElement_Standard
from objects.CSCG._3d.fields.vector.do.reconstruct.trace_element.trace_element_wise import OnTraceElement_TraceElementWise
from objects.CSCG._3d.fields.vector.do.reconstruct.trace_element.boundary_wise import OnTraceElement_BoundaryWise
from objects.CSCG._3d.fields.vector.do.reconstruct.trace_element.standard import OnTraceElement_Standard

class _3dCSCG_Vector_Do_Reconstruct(FrozenOnly):
    def __init__(self, vf):
        self._vf_ = vf
        self._on_mesh_element___for_standard_ = OnMeshElement_Standard(vf)
        self._on_trace_element___for_standard_ = OnTraceElement_Standard(vf)
        self._on_trace_element___for_trace_element_wise_ = OnTraceElement_TraceElementWise(vf)
        self._on_trace_element___for_boundary_wise_ = OnTraceElement_BoundaryWise(vf)
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
        ftype = self._vf_.ftype

        #------- deal with time --------------------------------------------------
        if time is None:
            pass
        else:
            self._vf_.current_time = time

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
            if ftype == 'standard':
                return self._on_trace_element___for_standard_(xi, eta, sigma, ravel, i)
            elif ftype == 'trace-element-wise':
                return self._on_trace_element___for_trace_element_wise_(xi, eta, sigma, ravel, i)
            elif ftype == 'boundary-wise':
                return self._on_trace_element___for_boundary_wise_(xi, eta, sigma, ravel, i)
            else:
                raise NotImplementedError()
        else:
            raise NotImplementedError()