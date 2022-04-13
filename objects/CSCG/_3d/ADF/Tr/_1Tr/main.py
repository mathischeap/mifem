
import sys
if './' not in sys.path: sys.path.append('./')

from objects.CSCG._3d.ADF.Tr.base.main import _3dCSCG_ADF_Tr_BASE


class _3dCSCG_1Tr_ADF(_3dCSCG_ADF_Tr_BASE):
    """
    Dual Tr 1-form.

    :param mesh:
    :param space:
    :param orientation:
    :param name:
    """
    def __init__(self, prime, mesh, space, orientation='outer', name=None):
        """"""
        if name is None: name = orientation + '-oriented-1-AD-Tr'
        super().__init__(3, mesh, space, prime, orientation, name)
        self._k_ = 1
        self.standard_properties.___PRIVATE_add_tag___('3dCSCG_standard_algebra_dual_1Tr')
        self.___PRIVATE_reset_cache___()
        self._freeze_self_()


    def ___PRIVATE_reset_cache___(self):
        """"""
        super().___PRIVATE_reset_cache___()






if __name__ == "__main__":
    # mpiexec -n 6 python objects\CSCG\_3d\ADF\Tr\_1Tr\main.py

    from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller#, ExactSolutionSelector

    # mesh = MeshGenerator('bridge_arch_cracked')([8,9,7], EDM='SWV0', show_info=True)
    mesh = MeshGenerator('crazy')([3, 3, 3], EDM=None, show_info=True)
    space = SpaceInvoker('polynomials')([2, 3, 2], show_info=True)
    # es = ExactSolutionSelector(mesh)('icpsNS:still', show_info=True)
    # ke0 = es.status.kinetic_energy(0)
    # he0 = es.status.helicity(0)

    FC = FormCaller(mesh, space)
    ad1 = FC('1-adTr', numbering_parameters={'scheme_name': 'Naive', })
