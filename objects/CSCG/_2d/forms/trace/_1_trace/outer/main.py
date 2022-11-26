# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
import sys
if './' not in sys.path: sys.path.append('/')
from objects.CSCG._2d.forms.trace._1_trace.base import _2dCSCG_1Trace
from objects.CSCG._2d.forms.trace._1_trace.outer.discretize.main import _2dCSCG_Outer0Trace_discretize


class _2dCSCG_1Trace_Outer(_2dCSCG_1Trace):
    """
    Outer trace 1-form.

    :param mesh:
    :param space:
    :param numbering_parameters:
    :param name:
    """
    def __init__(self, mesh, space, numbering_parameters='Naive', name='outer-oriented-1-trace-form'):
        super().__init__(mesh, space, True, 'outer', numbering_parameters, name)
        self._k_ = 1
        self.standard_properties.___PRIVATE_add_tag___('2dCSCG_trace_1form')
        self.RESET_cache()
        self._discretize_ = _2dCSCG_Outer0Trace_discretize(self)
        self._freeze_self_()

    def RESET_cache(self):
        super().RESET_cache()

    def ___PRIVATE_TW_FUNC_body_checker___(self, func_body):
        if func_body.__class__.__name__ == '_2dCSCG_ScalarField':
            assert func_body.mesh.domain == self.mesh.domain
            assert func_body.ndim == self.ndim == 2
        else:
            raise NotImplementedError()

    def discretize(self):
        return self._discretize_


    def reconstruct(self, xi, eta, ravel=False, key=None):
        """
        Do the reconstruction.

        :param xi: A 1d iterable object of floats between -1 and 1.
        :param eta: A 1d iterable object of floats between -1 and 1.
        :param bool ravel: (`default`:``False``) If we return 1d data?
        :param key: (`default`:``None``) Do the reconstruction for trace element named ``key``. if it is ``None``,
            then do it for all trace elements.
        :type xi: list, tuple, numpy.ndarray
        :type eta: list, tuple, numpy.ndarray
        :type key: str, None
        """
        raise NotImplementedError




if __name__ == '__main__':
    # mpiexec python _2dCSCG\form\trace\_1_trace_outer.py

    from objects.CSCG._2d.master import *
    # from mifem import read, save

    mesh = MeshGenerator('cic')([2,2])
    # mesh = MeshGenerator('crazy',)([2,2])
    space = SpaceInvoker('polynomials')([('Lobatto',2), ('Lobatto',3)])
    FC = FormCaller(mesh, space)

    # mesh.trace.visualize()

    t1o = FC('1-t-o')

    t1o.numbering.visualize(show_element_numbering=False, show_boundary_names=False)

    print(t1o.IS_hybrid)