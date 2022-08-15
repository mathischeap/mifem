# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/7/23 11:19
"""
import sys

if './' not in sys.path: sys.path.append('./')
from screws.freeze.main import FrozenOnly

import numpy as np
from tools.linear_algebra.nonlinear_system.solve.main import NonLinearSystem_Solve
from tools.linear_algebra.nonlinear_system.nonliner_terms import nLS_NonlinearTerms
from tools.linear_algebra.nonlinear_system.num import nLS_num
from tools.linear_algebra.nonlinear_system.do import nLS_DO
from tools.linear_algebra.nonlinear_system.customize import nLS_Customize


class NonLinearSystem(FrozenOnly):
    """"""
    def __init__(self, test_variables, linear_blocks, nonlinear_terms, unknown_variables, b):
        """

        Parameters
        ----------
        test_variables
        linear_blocks
        nonlinear_terms
        unknown_variables
        b
        """
        self.___PRIVATE_parse_test_and_unknown_variables___(test_variables, unknown_variables)

        self.___PRIVATE_parse_linear_blocks___(linear_blocks)
        self.___PRIVATE_parse_nonlinear_terms___(nonlinear_terms)
        self.___PRIVATE_parse_b___(b)
        assert self.regularity is not None, f"We must have assign it a regularity."

        self._solve_ = NonLinearSystem_Solve(self)
        self._num_ = nLS_num(self)
        self._do_ = nLS_DO(self)
        self._customize_ = nLS_Customize(self)
        self._freeze_self_()

    def __repr__(self):
        """"""
        return f"{self.shape}NonlinearSystem:{id(self)}"

    def ___PRIVATE_parse_test_and_unknown_variables___(self, test_variables, unknown_variables):
        """"""
        assert isinstance(test_variables, (tuple, list)), f"pls put test_variables in a list."
        assert isinstance(unknown_variables, (tuple, list)), f"pls put unknown variables in a list."

        if isinstance(test_variables, tuple): test_variables = list(test_variables)
        if isinstance(unknown_variables, tuple): unknown_variables = list(unknown_variables)

        for tv in test_variables:
            assert tv not in unknown_variables, f"test variables must be different from the unknowns."

        self._test_variables_ = test_variables
        self._unknown_variables_ = unknown_variables

        self._num_equations_ = 0
        for tv in test_variables:
            self._num_equations_ += tv.num.GLOBAL_dofs

        self._num_dofs_ = 0
        for uv in unknown_variables:
            self._num_dofs_ += uv.num.GLOBAL_dofs

    @property
    def unknown_variables(self):
        return self._unknown_variables_

    @property
    def test_variables(self):
        return self._test_variables_

    def ___PRIVATE_parse_linear_blocks___(self, A):
        """

        Parameters
        ----------
        A

        Returns
        -------

        """
        assert isinstance(A, (list, tuple)), f"Please put linear blocks in a list or a tuple."

        #---------- parse blocks ------------------------------------------------------
        shape0 = len(A)
        shape1 = None
        for i, Ai_ in enumerate(A):

            if isinstance(Ai_, tuple):
                Ai_ = list(Ai_)
                A[i] = Ai_

            if shape1 is None:
                shape1 = len(Ai_)
            else:
                assert shape1 == len(Ai_)

            for j, Aij in enumerate(Ai_):
                if Aij.__class__.__name__  == 'EWC_SparseMatrix':
                    pass
                elif Aij is None:
                    pass
                elif Aij == 0:
                    A[i][j] = None
                else:
                    raise Exception(f'A[{i}][{j}] = {Aij} is wrong.')

        self._shape_ = (shape0, shape1)
        assert (len(self.test_variables), len(self._unknown_variables_)) == self.shape, f"shape wrong."
        self._A_ = tuple(A)

    @property
    def shape(self):
        return self._shape_

    def ___PRIVATE_parse_nonlinear_terms___(self, terms):
        """"""
        assert isinstance(terms, (list, tuple)), f"pls put nonlinear terms in a list."
        assert len(terms) == self.shape[0], f"nonlinear terms length wrong."
        self.___nonlinear_terms___ = terms
        self._nonlinear_terms_ = nLS_NonlinearTerms(self)

        #----- parse regularity lhs: inhomogeneous_MDM -------------------------------------
        regularity = 'inhomogeneous_MDM'
        for i, row_terms in enumerate(terms):
            if row_terms == 0 or row_terms is None:
                pass
            elif row_terms.__class__.__name__ == 'MultiDimMatrix' and row_terms.IS.inhomogeneous:
                pass
            elif isinstance(row_terms, (list, tuple)):
                for term in row_terms:
                    if term.__class__.__name__ == 'MultiDimMatrix' and term.IS.inhomogeneous:
                        pass
                    else:
                        regularity = None
            else:
                regularity = None

        if regularity is not None:
            self._regularity_nonlinear_terms_ = regularity

            for tm_ind in self.nonlinear_terms:
                tm = self.nonlinear_terms[tm_ind]
                if tm not in (0, None):
                    assert tm.__class__.__name__ == 'MultiDimMatrix'

                    i = tm_ind[0]
                    tv = self.test_variables[i]
                    assert tv in tm.correspondence
                    for cor in tm.correspondence:
                        if cor is not tv:
                            assert cor in self.unknown_variables

            return

        #---------- new regularity check -----------------------------------------------

        #===============================================================================
        raise Exception() # find no suitable regularity, raise Exception.

    def ___PRIVATE_parse_b___(self, b):
        """

        Parameters
        ----------
        b

        Returns
        -------

        """
        assert isinstance(b, (list, tuple)), f"pls put right-hand-side vector in a list."
        if isinstance(b, tuple): b = list(b)
        #----------- save b -------------------------------------------------------------
        self._b_ = b

        #---- check regularity 1 of b --------------------------------------------------------
        regularity = 'EWC_cv'
        for bi in b:
            if bi.__class__.__name__  == 'EWC_ColumnVector':
                pass
            elif bi is None or bi == 0:
                pass
            else:
                regularity = None

        if regularity is not None:
            self._regularity_b_ = 'EWC_cv'
            return
        #------ check other regularities -----------------------------------------------------

        #=====================================================================================
        raise Exception() # find no suitable regularity for b, raise Error.

    @property
    def A(self):
        return self._A_

    @property
    def nonlinear_terms(self):
        return self._nonlinear_terms_

    @property
    def b(self):
        return self._b_

    @property
    def regularity(self):
        """"""
        return self._regularity_nonlinear_terms_, self._regularity_b_

    @property
    def solve(self):
        return self._solve_

    @property
    def num(self):
        return self._num_

    @property
    def do(self):
        return self._do_

    @property
    def customize(self):
        return self._customize_





if __name__ == '__main__':
    # mpiexec -n 4 python tools/linear_algebra/nonlinear_system/main.py
    from __init__ import cscg3

    mesh = cscg3.mesh('cuboid', region_layout=[2,2,2])([1,1,1])
    space = cscg3.space('polynomials')((2,2,2))
    FC = cscg3.form(mesh, space)

    def u(t,x,y,z): return np.sin(np.pi*x)*np.cos(2*np.pi*y)*np.cos(np.pi*z) + t
    def v(t,x,y,z): return np.cos(np.pi*x)*np.sin(np.pi*y)*np.cos(2*np.pi*z) + t
    def w(t,x,y,z): return np.cos(np.pi*x)*np.cos(np.pi*y)*np.sin(2*np.pi*z) + t

    velocity = FC('vector', (u,v,w))
    U = FC('scalar', u)
    V = FC('scalar', v)
    W = FC('scalar', w)

    w = FC('1-f', is_hybrid=False)
    o = FC('1-f', is_hybrid=False)
    u = FC('2-f', is_hybrid=False)
    v = FC('2-f', is_hybrid=False)

    MDM = w.special.cross_product_2f__ip_2f(u, v, output='MDM')

    M2 = u.matrices.mass
    M1 = w.matrices.mass

    w.TW.func.do.set_func_body_as(velocity)
    w.TW.current_time = 0
    w.TW.do.push_all_to_instant()
    w.discretize()

    u.TW.func.do.set_func_body_as(velocity)
    u.TW.current_time = 0
    u.TW.do.push_all_to_instant()
    u.discretize()

    b0 = u.cochain.EWC
    b1 = w.cochain.EWC

    nLS = NonLinearSystem([v, o], ([M2, 0 ],
                                   [ 0, M1]), [MDM, 0],
                          [u,w], [b0, b1])
    print(nLS.nonlinear_terms.num)