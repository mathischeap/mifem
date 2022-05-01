
from screws.freeze.base import FrozenOnly
from objects.CSCG._3d.forms.standard.base.dofs.dof.do.generate_plot_data._0sf import GPD_0SF
from objects.CSCG._3d.forms.standard.base.dofs.dof.do.generate_plot_data._1sf import GPD_1SF
from objects.CSCG._3d.forms.standard.base.dofs.dof.do.generate_plot_data._2sf import GPD_2SF
from objects.CSCG._3d.forms.standard.base.dofs.dof.do.generate_plot_data._3sf import GPD_3SF

class _3dCSCG_SF_dof_DO(FrozenOnly):
    """"""
    def __init__(self, dof):
        """"""
        self._dof_ = dof
        self._freeze_self_()

    def generate_plot_data_from_element(self, mesh_element, *args, **kwargs):
        """

        Parameters
        ----------
        mesh_element : int
            We make plotting data considering the dof is in this mesh-element.
        args :
        kwargs :

        Returns
        -------

        """
        if self._dof_._sf_.k == 0:
            return GPD_0SF(self._dof_)(mesh_element, *args, **kwargs)
        elif self._dof_._sf_.k == 1:
            return GPD_1SF(self._dof_)(mesh_element, *args, **kwargs)
        elif self._dof_._sf_.k == 2:
            return GPD_2SF(self._dof_)(mesh_element, *args, **kwargs)
        elif self._dof_._sf_.k == 3:
            return GPD_3SF(self._dof_)(mesh_element, *args, **kwargs)
        else:
            raise Exception()