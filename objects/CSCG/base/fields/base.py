# -*- coding: utf-8 -*-
from root.config.main import *
from objects.base.fields.base import FiledBase


class CSCG_Continuous_FORM_BASE(FiledBase):
    """"""
    def __init_subclass__(cls, ndim=None):
        cls.___ndim___ = ndim

    def __init__(self, mesh, ftype, valid_time):
        super(CSCG_Continuous_FORM_BASE, self).__init__(mesh, valid_time)
        self._visualize_ = None
        self._CMB_ = None
        self.standard_properties.___PRIVATE_add_tag___('CSCG_field')

        assert ftype in ("standard", "boundary-wise", "trace-element-wise"), \
            f"ftype={ftype} wrong." # allowed ftype.


    @property
    def ndim(self):
        return self.___ndim___

    @property
    def mesh(self):
        """
        The mesh.

        The mesh will only be cleared when we save a form have has this vector field as property.
        """
        return self._mesh_

    @mesh.setter
    def mesh(self, mesh):
        assert mesh.__class__.__name__ == f'_{self.ndim}dCSCG_Mesh', f"Need a _{self.ndim}dCSCG_Mesh."
        self._mesh_ = mesh

    # --------------------------------------------------------------------------------------
    @property
    def visualize(self):
        return self._visualize_

    @property
    def covered_mesh_boundaries(self):
        """The covered mesh boundaries: This continuous field is valid (defined)
        on these mesh boundaries (NOT domain boundaries). Same in all cores!
        So even some mesh boundaries is not included on the local mesh elements,
        we still return them with this property."""

        if self._CMB_ is not None: return  self._CMB_

        AMBs = self.mesh.boundaries.names # All Mesh Boundaries

        if self.ftype == 'standard':
            self._CMB_ = AMBs

        elif self.ftype == 'boundary-wise':
            self._CMB_ = list()
            for bn in self.func:
                assert bn in AMBs
                self._CMB_.append(bn)
        elif self.ftype == 'trace-element-wise':
            _CMB_ = list()
            for i in self.func: # valid local trace elements
                te = self.mesh.trace.elements[i]
                if te.whether.on_mesh_boundary:
                    omb = te.on_mesh_boundary
                    if omb not in _CMB_:
                        _CMB_.append(omb)
                else:
                    pass

            _CMB_ = COMM.gather(_CMB_, root=MASTER_RANK)
            if RANK == MASTER_RANK:
                self._CMB_ = set()
                for __ in _CMB_:
                    self._CMB_.update(__)
                self._CMB_ = list(self._CMB_)
            else:
                self._CMB_ = None

            self._CMB_ = COMM.bcast(self._CMB_, root=MASTER_RANK)

        else:
            raise NotImplementedError(f"Can not deal with {self.ftype} ftype.")

        return self._CMB_

    #====================== BELOW: children must have these methods =======================
    def ___PRIVATE_set_func___(self, *args, **kwargs):
        raise NotImplementedError()

    def ___DO_evaluate_func_at_time___(self, *args, **kwargs):
        raise NotImplementedError()

    @property
    def shape(self):
        raise NotImplementedError()

    def reconstruct(self, *args, **kwargs):
        raise NotImplementedError()

    @property
    def do(self):
        raise NotImplementedError()