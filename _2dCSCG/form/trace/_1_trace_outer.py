# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
import sys
if './' not in sys.path: sys.path.append('./')
from _2dCSCG.form.trace.main import _2dCSCG_Standard_Trace


class _1Trace_Outer(_2dCSCG_Standard_Trace):
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
        self._freeze_self_()

    def RESET_cache(self):
        self.___cache_DISCRETIZE_STANDARD___ = None
        super().RESET_cache()

    def ___TW_FUNC_body_checker___(self, func_body):
        if func_body.__class__.__name__ == '_2dCSCG_ScalarField':
            assert func_body.mesh.domain == self.mesh.domain
            assert func_body.ndim == self.ndim == 2
        else:
            raise NotImplementedError()

    def discretize(self, update_cochain=True, **kwargs):
        """
        Do the discretization.

        :param bool update_cochain: Whether we update the cochain if the trace form.
        :param kwargs: Keywords arguments to be passed to particular discretization schemes.
        :return: The cochain corresponding to the particular discretization scheme.
        """
        if self.func.ftype == 'standard':
            return self.___PRIVATE_discretize_standard_ftype___(
                update_cochain=update_cochain, **kwargs)
        else:
            raise NotImplementedError()

    def ___PRIVATE_discretize_standard_ftype___(self, update_cochain=True, quad_degree=None):
        """We will discretize the a standard scalar field to all trace elements."""
        raise NotImplementedError()

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

    def ___DO_resemble___(self, obj_or_filename):
        """

        :param obj_or_filename:
        :return:
        """
        raise NotImplementedError()




if __name__ == '__main__':
    # mpiexec python _2dCSCG\form\trace\_1_trace_outer.py

    from _2dCSCG.main import *
    # from mifem import read, save

    mesh = MeshGenerator('cic')([2,2])
    # mesh = MeshGenerator('crazy',)([2,2])
    space = SpaceInvoker('polynomials')([('Lobatto',2), ('Lobatto',3)])
    FC = FormCaller(mesh, space)

    # mesh.trace.visualize()

    t1o = FC('1-t-o')

    t1o.numbering.visualize(show_element_numbering=False, show_boundary_names=False)

    print(t1o.IS_hybrid)