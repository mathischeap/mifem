


from screws.freeze.base import FrozenOnly


from objects.CSCG._2d.fields.vector.do.reconstruct.mesh_element.standard import OnMeshElement_for_Standard




class _2dCSCG_Vector_Do_Reconstruct(FrozenOnly):
    """"""
    def __init__(self, vf):
        self._vf_ = vf
        self._on_mesh_element___for_standard_ = OnMeshElement_for_Standard(vf)
        self._freeze_self_()

    def __call__(self, xi, eta, time=None, ravel=False, i=None, where=None):
        """

        :param xi:
        :param eta:
        :param time:
        :param ravel:
        :param i:
        :param where:
        :return:
        """
        ftype = self._vf_.ftype

        # --- deal with time ---------------------------------------------------------
        if time is None:
            pass
        else:
            self._vf_.current_time = time

        # parse where when it is None --------------------------------------------
        if where is None:
            if ftype == 'standard':
                where = 'mesh-element'
            else:
                raise NotImplementedError(f"please set default `where` for {ftype} 2dCSCG scalar field.")
        else:
            pass

        #--------------------------------------------------------------------------
        if where == 'mesh-element':
            if ftype == 'standard':
                return self._on_mesh_element___for_standard_(xi, eta, ravel, i)
            else:
                raise NotImplementedError()
        else:
            raise NotImplementedError()