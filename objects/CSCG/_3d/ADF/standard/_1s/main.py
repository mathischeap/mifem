import sys
if './' not in sys.path: sys.path.append('./')


from objects.CSCG._3d.ADF.standard.base.main import _3dCSCG_Algebra_DUAL_Standard_Form


class _3dCSCG_S1_ADF(_3dCSCG_Algebra_DUAL_Standard_Form):
    """
    Standard a dual 1-form.

    :param mesh:
    :param space:
    :param orientation:
    :param name:
    """
    def __init__(self, prime, mesh, space, orientation='outer', name=None):
        if name is None: name = orientation + '-oriented-1-ADF'
        super().__init__(3, mesh, space, prime, orientation, name)
        self._k_ = 1
        self.standard_properties.___PRIVATE_add_tag___('3dCSCG_standard_algebra_dual_1form')
        self.___PRIVATE_reset_cache___()
        self._freeze_self_()




    def ___PRIVATE_reset_cache___(self):
        """"""
        super().___PRIVATE_reset_cache___()




if __name__ == "__main__":
    # mpiexec -n 6 python _3dCSCG\ADF\standard\_1_AD_form.py
    from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller
    import numpy as np

    # mesh = MeshGenerator('bridge_arch_cracked')([8,9,7], EDM='SWV0', show_info=True)
    mesh = MeshGenerator('crazy', c = 0.)([10, 10, 10], EDM=None, show_info=True)
    space = SpaceInvoker('polynomials')([('Lobatto',3), ('Lobatto',3), ('Lobatto',3)], show_info=True)


    def u(t, x, y, z):
        return np.cos(2*np.pi*x) * np.sin(2*np.pi*y) * np.sin(2*np.pi*z) + t
    def v(t, x, y, z):
        return np.sin(2*np.pi*x) * np.cos(2*np.pi*y) * np.sin(2*np.pi*z) + t
    def w(t, x, y, z):
        return np.sin(2*np.pi*x) * np.sin(2*np.pi*y) * np.cos(2*np.pi*z) + t

    def p(t, x, y, z):
        return - 6 * np.pi * np.sin(2*np.pi*x) * np.sin(2*np.pi*y) * np.sin(2*np.pi*z) + 0 * t

    FC = FormCaller(mesh, space)
    df1 = FC('1-adf', numbering_parameters={'scheme_name': 'Naive', })
    dt0 = FC('0-adt', numbering_parameters={'scheme_name': 'Naive', })

    scalar = FC('scalar', p)
    vector = FC('vector', (u, v, w))

    df1.prime.TW.func.do.set_func_body_as(vector)
    df1.prime.TW.current_time = 0
    df1.prime.TW.do.push_all_to_instant()
    df1.prime.do.discretize()

    dt0.prime.TW.func.do.set_func_body_as(vector)
    dt0.prime.TW.current_time = 0
    dt0.prime.TW.do.push_all_to_instant()
    dt0.prime.do.discretize()

    df0 = df1.coboundary(dt0)

    # df0 = FC('0-adf')
    df0.prime.TW.func.do.set_func_body_as(scalar)
    df0.prime.TW.current_time = 0
    df0.prime.TW.do.push_all_to_instant()
    # df0.prime.do.discretize()

    # df0.prime.visualize()

    print(df1.prime.error.L())
    print(df0.prime.error.L())