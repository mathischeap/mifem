import sys
if './' not in sys.path: sys.path.append('./')


from _3dCSCG.ADF.standard.main import _3dCSCG_Algebra_DUAL_Standard_Form


class _2_Algebra_DUAL_Form(_3dCSCG_Algebra_DUAL_Standard_Form):
    """
    Standard a dual 2-form.

    :param mesh:
    :param space:
    :param orientation:
    :param name:
    """
    def __init__(self, prime, mesh, space, orientation='outer', name=None):
        if name is None: name = orientation + '-oriented-2-ADF'
        super().__init__(3, mesh, space, orientation, name)
        self._prime_ = prime
        self._k_ = 2
        self.standard_properties.___PRIVATE_add_tag___('3dCSCG_standard_algebra_dual_2form')
        self.RESET_cache()
        self._freeze_self_()




    def RESET_cache(self):
        """"""
        super().RESET_cache()




if __name__ == "__main__":
    # mpiexec -n 6 python _3dCSCG\ADF\standard\_2_AD_form.py
    from _3dCSCG.main import MeshGenerator, SpaceInvoker, FormCaller
    import numpy as np

    # mesh = MeshGenerator('bridge_arch_cracked')([8,9,7], EDM='SWV0', show_info=True)
    mesh = MeshGenerator('crazy')([10, 10, 10], EDM=None, show_info=True)
    space = SpaceInvoker('polynomials')([('Lobatto',3), ('Lobatto',3), ('Lobatto',3)], show_info=True)



    def F(t, x, y, z):
        return np.sin(2*np.pi*x) * np.sin(2*np.pi*y) * np.sin(2*np.pi*z) + t
    def G(t, x, y, z):
        return np.sin(2*np.pi*x) * np.cos(2*np.pi*y) * np.sin(2*np.pi*z) + t
    def H(t, x, y, z):
        return np.sin(2*np.pi*x) * np.sin(2*np.pi*y) * np.cos(np.pi*z) + t

    def Fy(t, x, y, z):
        return 2 * np.pi * np.sin(2*np.pi*x) * np.cos(2*np.pi*y) * np.sin(2*np.pi*z) + 0*t
    def Fz(t, x, y, z):
        return 2 * np.pi * np.sin(2*np.pi*x) * np.sin(2*np.pi*y) * np.cos(2*np.pi*z) + 0*t

    def Gz(t, x, y, z):
        return 2 * np.pi * np.sin(2*np.pi*x) * np.cos(2*np.pi*y) * np.cos(2*np.pi*z) + 0*t
    def Gx(t, x, y, z):
        return 2 * np.pi * np.cos(2*np.pi*x) * np.cos(2*np.pi*y) * np.sin(2*np.pi*z) + 0*t

    def Hx(t, x, y, z):
        return 2 * np.pi * np.cos(2*np.pi*x) * np.sin(2*np.pi*y) * np.cos(np.pi*z) + 0*t
    def Hy(t, x, y, z):
        return 2 * np.pi * np.sin(2*np.pi*x) * np.cos(2*np.pi*y) * np.cos(np.pi*z) + 0*t


    def Cx(t, x, y, z): return Hy(t, x, y, z) - Gz(t, x, y, z)
    def Cy(t, x, y, z): return Fz(t, x, y, z) - Hx(t, x, y, z)
    def Cz(t, x, y, z): return Gx(t, x, y, z) - Fy(t, x, y, z)



    FC = FormCaller(mesh, space)
    df2 = FC('2-adf', numbering_parameters={'scheme_name': 'Naive', })
    dt1 = FC('1-adt', numbering_parameters={'scheme_name': 'Naive', })

    V = FC('vector', (F, G, H))
    C = FC('vector', (Cx, Cy, Cz))

    df2.prime.TW.func.DO.set_func_body_as(V)
    df2.prime.TW.current_time = 0
    df2.prime.TW.DO.push_all_to_instant()
    df2.prime.DO.discretize()

    dt1.prime.TW.func.DO.set_func_body_as(C)
    dt1.prime.TW.current_time = 0
    dt1.prime.TW.DO.push_all_to_instant()
    dt1.prime.DO.discretize()


    # df1 = FC('1-adf')
    df1 = df2.coboundary(dt1)

    df1.prime.TW.func.DO.set_func_body_as(C)
    df1.prime.TW.current_time = 0
    df1.prime.TW.DO.push_all_to_instant()
    # df1.prime.DO.discretize()


    print(df2.prime.error.L())
    print(df1.prime.error.L())