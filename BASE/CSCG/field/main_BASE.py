

from SCREWS.frozen import FrozenClass

from root.config import *

class CSCG_Continuous_FORM_BASE(FrozenClass):
    """"""
    def __init_subclass__(cls, ndim=None):
        cls.___ndim___ = ndim

    def __init__(self, mesh, ftype, valid_time):
        self._mesh_ = mesh
        self._func_ = None
        self._ftype_ = None
        self._current_time_ = None
        self._visualize_ = None
        self.standard_properties.___PRIVATE_add_tag___('CSCG_field')
        self._CMB_ = None
        self.___PRIVATE_set_vt_to___(valid_time)

        assert ftype in ("standard", "boundary-wise", "trace-element-wise"), f"ftype={ftype} wrong."

    # ... most standard properties ...
    @property
    def func(self):
        """The function body: same in all cores!"""
        return self._func_

    @property
    def ftype(self):
        """The function type: same in all cores!"""
        return self._ftype_


    # ....
    @property
    def ndim(self):
        return self.___ndim___

    @property
    def shape(self):
        raise NotImplementedError()

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




    # ------------------------------------------------------------------
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
                if te.IS_on_mesh_boundary:
                    omb = te.on_mesh_boundary
                    if omb not in _CMB_:
                        _CMB_.append(omb)
                else:
                    pass

            _CMB_ = cOmm.gather(_CMB_, root=mAster_rank)
            if rAnk == mAster_rank:
                self._CMB_ = set()
                for __ in _CMB_:
                    self._CMB_.update(__)
                self._CMB_ = list(self._CMB_)
            else:
                self._CMB_ = None

            self._CMB_ = cOmm.bcast(self._CMB_, root=mAster_rank)

        else:
            raise NotImplementedError(f"Can not deal with {self.ftype} ftype.")

        return self._CMB_




    # ...........CURRENT TIME which is related to valid time ..........................................
    @property
    def current_time(self):
        """The current time. If push to instant, we use this time."""
        assert self._current_time_ is not None, "current_time is None, set it firstly."
        return self._current_time_

    @current_time.setter
    def current_time(self, current_time):
        assert isinstance(current_time, (float,int))
        self.___PRIVATE_check_ct_in_vt___(current_time)
        self._current_time_ = current_time


    #----------- VALID TIME -----------------------------------------------------------------
    @property
    def valid_time(self):
        """``current time`` must be with in ``valid_time``.

        Some continuous forms are only valid for some certain times.

        `valid_time` cannot be changed once we have define the instance.
        """
        return self._valid_time_

    def ___PRIVATE_set_vt_to___(self, valid_time):
        """
        :param valid_time:
            None: It can be everything and be changed whenever you want.
            'valid_only_at_its_first_instant': as it says...
        :return:
        """
        if valid_time is None:
            pass
        elif valid_time == 'valid_only_at_its_first_instant':
            # as the string says, we can only set ``current_time`` once.
            pass
        else:
            raise Exception(f'valid_time = {valid_time} format wrong.')
        self._valid_time_ = valid_time

    def ___PRIVATE_check_ct_in_vt___(self, ct):
        if self.valid_time is None:
            # valid at any given time instant (current time).
            return
        elif self.valid_time == 'valid_only_at_its_first_instant':
            # current time is only valid at its first set instant.
            if self._current_time_ is None:
                pass
            elif ct == self._current_time_:
                pass
            else:
                raise Exception(f"time is `valid_only_at_its_first_instant`. "
                                f"current_time = {self._current_time_}, "
                                f"can not change it to {ct}.")
        else:
            raise Exception(f'current_time = {ct} can not be checked')

    #=======================================================================================




    def reconstruct(self, *args, **kwargs):
        raise NotImplementedError()

    def ___DO_evaluate_func_at_time___(self, *args, **kwargs):
        raise NotImplementedError()

    def ___PRIVATE_set_func___(self, *args, **kwargs):
        raise NotImplementedError()

    @property
    def DO(self):
        raise NotImplementedError()
