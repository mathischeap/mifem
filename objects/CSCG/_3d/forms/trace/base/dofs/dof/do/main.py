


from screws.freeze.base import FrozenOnly


from objects.CSCG._3d.forms.trace.base.dofs.dof.do.generate_plot_data._0where import _3dCSCG_T0F_DOF_Where
from objects.CSCG._3d.forms.trace.base.dofs.dof.do.generate_plot_data._1where import _3dCSCG_T1F_DOF_Where
from objects.CSCG._3d.forms.trace.base.dofs.dof.do.generate_plot_data._2where import _3dCSCG_T2F_DOF_Where

class _3dCSCG_TF_dof_DO(FrozenOnly):
    """"""
    def __init__(self, dof):
        """"""
        self._dof_ = dof
        self._freeze_self_()

    def generate_plot_data(self, *args, **kwargs):
        """

        Parameters
        ----------
        args :
        kwargs :

        Returns
        -------

        """
        k = self._dof_._tf_.k
        if k == 0:
            return _3dCSCG_T0F_DOF_Where(self._dof_)(*args, **kwargs)
        elif k == 1:
            return _3dCSCG_T1F_DOF_Where(self._dof_)(*args, **kwargs)
        elif k == 2:
            return _3dCSCG_T2F_DOF_Where(self._dof_)(*args, **kwargs)
        else:
            raise Exception()