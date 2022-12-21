# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 11/28/2022 2:29 PM
"""
from components.freeze.main import FrozenOnly


class _3dCSCG_LocalTrace_Do(FrozenOnly):
    """"""

    def __init__(self, ltf):
        """"""
        self._ltf_ = ltf
        self._freeze_self_()

    def p_refine(self, p=(1, 1, 1)):
        """Return another instance whose degree is higher of p.

        For example, ltf.p = (2, 3, 4), ltf.do.p_refine(p=(3, 1, 2) gives another instance
        of the class of ltf, and its space degree is (2+3, 3+1, 4+2).

        For the time being, we do not update the cochain of the new form. So its cochain.local
        is None.

        Parameters
        ----------
        p

        Returns
        -------

        """
        assert len(p) == 3 and all([isinstance(_, int) for _ in p]), \
            f"p = {p} is wring, p must be a tuple or list of 3 integers (no need to be positive)."

        new_space = self._ltf_.space.do.refine(p=p)
        numbering_parameters = {'scheme_name': self._ltf_.numbering._scheme_name_}
        numbering_parameters.update(self._ltf_.numbering._parameters_)

        return self._ltf_.__class__(
            self._ltf_.mesh, new_space,
            orientation=self._ltf_._orientation_,
            numbering_parameters=numbering_parameters,
            name='p-refined-' + self._ltf_._orientation_ + '-0-ltf'
        )

    def evaluate_basis_at_meshgrid(self, xi, eta, sigma, compute_xietasigma=True):
        """
        Evaluate the basis functions on ``meshgrid(xi, eta, sigma)``.

        :param xi: A 1d iterable object of floats between -1 and 1.
        :param eta: A 1d iterable object of floats between -1 and 1.
        :param sigma: A 1d iterable object of floats between -1 and 1.
        :type xi: list, tuple, numpy.ndarray
        :type eta: list, tuple, numpy.ndarray
        :type sigma: list, tuple, numpy.ndarray
        :param bool compute_xietasigma: (`default`:``True``) If we compute the
            ``meshgrid(xi, eta, sigma, indexing='ij')``.
        :returns: A tuple of outputs:

            1. (None, tuple) -- ``(xi, eta, sigma)`` after ``meshgrid`` and ``ravel('F')``.
            2. tuple -- The evaluated basis functions.
        """
        return self._ltf_.space.do.evaluate_local_trace_basis_at_meshgrid(
            self._ltf_.k, xi, eta, sigma, compute_xietasigma=compute_xietasigma)

    def make_reconstruction_matrix_on_grid(self, xi, eta, sigma, element_range=None):
        """Make the reconstruction matrices for all mesh elements.

        These matrices are stored in a dict whose keys are the numbers of mesh elements and values
        are the local reconstruction matrices.

        Let `RM` be the reconstruction matrix (or the tuple of three matrices).
        If we want to do the local reconstruction, we do

            RM[i] @ f.cochain.local[i]

        and we will get the reconstructions of the form `f` on `meshgrid(xi, eta, sigma)` in mesh-element
        #i. And if `f` is a scalar form, we get a 1d array. And if `f` is a vector form, we get a
        tuple of three 1d arrays (its three components along x, y, z directions.)

        :param xi: 1d array
        :param eta: 1d array
        :param sigma: 1d array
        :param element_range:
            We are going to construct matrices for these mesh elements. It can be one of
                1) None: for all local elements
                2) 'mesh boundary': those local elements are attached to mesh boundary.

        """
        return self._ltf_.reconstruct.___PrLT_make_reconstruction_matrix_on_grid___(
            xi, eta, sigma, element_range=element_range
        )
