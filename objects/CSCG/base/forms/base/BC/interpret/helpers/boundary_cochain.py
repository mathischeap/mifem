# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 12/15/2022 4:42 PM
"""
from scipy.sparse import csc_matrix, csr_matrix
from components.freeze.main import FrozenOnly
import numpy as np


class CSCG_FORM_BC_Interpret_BoundaryCochain(FrozenOnly):
    """It will follow the s.BC.CF and s.BC.boundaries in real time."""

    def __init__(self, f):
        """"""
        self._f_ = f
        self._mesh_ = f.mesh
        self._elements_ = f.mesh.elements
        self.___empty___ = csc_matrix((f.num.basis, 1))
        self.___boundary_cochain_cache___ = None
        self._representing_CF_ = -1  # will store id, so use int
        self._representing_ct_ = None
        self._representing_boundaries_ = -1  # will store id, so use int
        self._freeze_self_()

    def ___Pr_make_mesh_element_wise_boundary_local_dofs___(self):
        """"""
        assert self._f_.BC.CF is not None, f"BC boundary is not empty, so please set boundary function, i.e. `BC.CF`."
        indicator, local_cochain = self._f_.discretize(target='BC')
        self.___boundary_cochain_cache___ = dict()  # renew the `boundary_cochain_cache`.

        if indicator == 'mesh-element-side-wise local cochain':
            # will only look at cochains for dofs on mesh element side (boundary of the mesh).
            iEP = self._f_.BC._involved_element_parts_

            for e_p in iEP:
                e, side = int(e_p[:-1]), e_p[-1]
                assert e in local_cochain, \
                    f"element {e} in not in the local cochain, most likely," \
                    f"the boundaries in the func do not cover `BC.boundaries`: " \
                    f"{self._f_.BC.boundaries}."

                if e not in self.___boundary_cochain_cache___:
                    self.___boundary_cochain_cache___[e] = np.zeros(self._f_.num.basis)
                else:
                    pass

                if self._mesh_.ndim == 3:
                    dofs = self._f_.numbering.do.\
                        find.local_dofs_on_element_side(side)
                elif self._mesh_.ndim == 2:
                    dofs = self._f_.numbering.do.\
                        find.local_dofs_on_element_edge(side)
                else:
                    raise Exception()

                self.___boundary_cochain_cache___[e][dofs] = local_cochain[e][side]

        else:
            raise NotImplementedError(f"indicator={indicator} is not implemented.")

        # we update representing information below.
        self._representing_CF_ = id(self._f_.BC.CF)
        # we use `CF.current_time` here as if we reach here but do not see current time, it should raise error.
        self._representing_ct_ = self._f_.BC.CF.current_time
        # we update representing boundaries as well
        self._representing_boundaries_ = id(self._f_.BC.boundaries)

    @property
    def ___Pr_BcCo___(self):
        if self.___boundary_cochain_cache___ is None:
            self.___Pr_make_mesh_element_wise_boundary_local_dofs___()
        else:
            # the following check is OK, since when CF is None, we will never call `___Pr_BcCo___`.
            if self._representing_CF_ == id(self._f_.BC.CF) and \
                    self._representing_ct_ == self._f_.BC.CF.current_time and \
                    self._representing_boundaries_ == id(self._f_.BC.boundaries):
                # we use `CF.current_time` here as if we reach here but do not see current time, it should raise error.
                pass  # as nothing has changed, so we use the cached instance.
            else:
                self.___Pr_make_mesh_element_wise_boundary_local_dofs___()

        return self.___boundary_cochain_cache___

    def __call__(self, i):
        """Return the boundary local cochains in real time (following the current `BC.CF` and `BC.boundaries`).
        for mesh-element i. For those dofs not locating on the mesh boundary, we set 0 for them.

        Parameters
        ----------
        i

        Returns
        -------

        """
        if self._f_.BC.boundaries == list() or \
           self._f_.BC.boundaries is None:

            return self.___empty___

        else:
            element = self._elements_[i]
            if element.whether.internal:
                # cause this is a realtime function, we always check self._f_.BC.
                return self.___empty___

            else:
                if i in self.___Pr_BcCo___:
                    return csr_matrix(self.___Pr_BcCo___[i]).T
                else:
                    return self.___empty___
