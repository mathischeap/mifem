
"""
For standard forms only.

"""

import sys
if './' not in sys.path: sys.path.append('./')
import os
from root.config.main import *
from scipy.sparse import linalg as spspalinalg
from _3dCSCG.main import MeshGenerator, SpaceInvoker, FormCaller, ExactSolutionSelector
from root.mifem.save import save
import random
from screws.exceptions import ThreeDimensionalTransfiniteInterpolationError
from _3dCSCG.tests.random_objects import random_3D_mesh_and_space_of_total_load_around
from _3dCSCG.tests.random_objects import random_3D_FormCaller_of_total_load_around


def test_Form_NO1_discretization_and_reconstruction():
    """"""
    if rAnk == mAster_rank:
        print("*** [test_Form_NO1_discretization_and_reconstruction] ...... ", flush=True)

    def u(t, x, y, z): return np.cos(np.pi*x) + np.sin(np.pi*y) * np.sin(np.pi*z-0.125)**2 + t/2
    def v(t, x, y, z): return np.sin(np.pi*x) + np.sin(np.pi*y) * np.sin(np.pi*z-0.125)**2 + t/2
    def w(t, x, y, z): return np.sin(np.pi*x) + np.cos(np.pi*y) * np.cos(np.pi*z-0.125)**2 + t
    def p(t, x, y, z): return np.cos(np.pi*x) + np.sin(np.pi*y) * np.sin(np.pi*z-0.125)**2 + t/2

    try:
        mesh = MeshGenerator('LDC', l=1, w=1.2, h=1.5)([f'Lobatto:{3}', f'Lobatto:{4}', f'Lobatto:{5}'], EDM='debug')
        space = SpaceInvoker('polynomials')([('Lobatto', 5), ('Lobatto', 5), ('Lobatto', 5)])
        FC = FormCaller(mesh, space)
        scalar = FC('scalar', p)
        vector = FC('vector', (u,v,w))
        f0 = FC('0-f', is_hybrid=False)
        f1 = FC('1-f', is_hybrid=False)
        f2 = FC('2-f', is_hybrid=False)
        f3 = FC('3-f', is_hybrid=False)
        f0.TW.func.body = scalar
        f0.TW.do.push_all_to_instant(0)
        f0.discretize()
        assert f0.error.L() < 0.00003
        f1.TW.func.body = vector
        f1.TW.do.push_all_to_instant(0)
        f1.discretize()
        assert f1.error.L() < 0.0004
        f2.TW.func.body = vector
        f2.TW.do.push_all_to_instant(0)
        f2.discretize()
        assert f2.error.L() < 0.0005
        f3.TW.func.body = scalar
        f3.TW.do.push_all_to_instant(0)
        f3.discretize()
        assert f3.error.L() < 0.0004
    except ThreeDimensionalTransfiniteInterpolationError:
        if rAnk == mAster_rank:
            print("    skip LDC mesh ...... ", flush=True)
        pass

    mesh = MeshGenerator('crazy_periodic', c=0.1)([3, 4, 5], EDM='debug')
    space = SpaceInvoker('polynomials')([('Lobatto', 5), ('Lobatto', 4), ('Lobatto', 3)])
    FC = FormCaller(mesh, space)

    scalar = FC('scalar', p)
    vector = FC('vector', (u,v,w))

    f0 = FC('0-f', is_hybrid=False)
    f1 = FC('1-f', is_hybrid=False)
    f2 = FC('2-f', is_hybrid=False)
    f3 = FC('3-f', is_hybrid=False)

    f0.TW.func.body = scalar
    f0.TW.do.push_all_to_instant(0)
    f0.discretize()
    assert f0.error.L() < 0.0004
    f1.TW.func.body = vector
    f1.TW.do.push_all_to_instant(0)
    f1.discretize()
    assert f1.error.L() < 0.004
    f2.TW.func.body = vector
    f2.TW.do.push_all_to_instant(0)
    f2.discretize()
    assert f2.error.L() < 0.0041
    f3.TW.func.body = scalar
    f3.TW.do.push_all_to_instant(0)
    f3.discretize()
    assert f3.error.L() < 0.003

    f0.cochain.globe = f0.cochain.globe
    f1.cochain.globe = f1.cochain.globe
    f2.cochain.globe = f2.cochain.globe
    f3.cochain.globe = f3.cochain.globe
    f0.TW.func.body = scalar
    f0.TW.do.push_all_to_instant(0)
    f0.discretize()
    assert f0.error.L() < 0.0004
    f1.TW.func.body = vector
    f1.TW.do.push_all_to_instant(0)
    f1.discretize()
    assert f1.error.L() < 0.004
    f2.TW.func.body = vector
    f2.TW.do.push_all_to_instant(0)
    f2.discretize()
    assert f2.error.L() < 0.0041
    f3.TW.func.body = scalar
    f3.TW.do.push_all_to_instant(0)
    f3.discretize()
    assert f3.error.L() < 0.003

    f0 = FC('0-f', is_hybrid=True)
    f1 = FC('1-f', is_hybrid=True)
    f2 = FC('2-f', is_hybrid=True)
    f3 = FC('3-f', is_hybrid=True)

    f0.TW.func.body = scalar
    f0.TW.do.push_all_to_instant(0)
    f0.discretize()
    assert f0.error.L() < 0.0004
    f1.TW.func.body = vector
    f1.TW.do.push_all_to_instant(0)
    f1.discretize()
    assert f1.error.L() < 0.004
    f2.TW.func.body = vector
    f2.TW.do.push_all_to_instant(0)
    f2.discretize()
    assert f2.error.L() < 0.0041
    f3.TW.func.body = scalar
    f3.TW.do.push_all_to_instant(0)
    f3.discretize()
    assert f3.error.L() < 0.003

    f0.cochain.globe = f0.cochain.globe
    f1.cochain.globe = f1.cochain.globe
    f2.cochain.globe = f2.cochain.globe
    f3.cochain.globe = f3.cochain.globe
    f0.TW.func.body = scalar
    f0.TW.do.push_all_to_instant(0)
    f0.discretize()
    assert f0.error.L() < 0.0004
    f1.TW.func.body = vector
    f1.TW.do.push_all_to_instant(0)
    f1.discretize()
    assert f1.error.L() < 0.004
    f2.TW.func.body = vector
    f2.TW.do.push_all_to_instant(0)
    f2.discretize()
    assert f2.error.L() < 0.0041
    f3.TW.func.body = scalar
    f3.TW.do.push_all_to_instant(0)
    f3.discretize()
    assert f3.error.L() < 0.003


    mesh = MeshGenerator('bridge_arch_cracked',)([3,2,4], EDM='debug')
    space = SpaceInvoker('polynomials')([('Lobatto', 3), ('Lobatto', 4), ('Lobatto', 2)])
    FC = FormCaller(mesh, space)

    scalar = FC('scalar', p)
    vector = FC('vector', (u,v,w))
    f0 = FC('0-f', is_hybrid=False)
    f1 = FC('1-f', is_hybrid=False)
    f2 = FC('2-f', is_hybrid=False)
    f3 = FC('3-f', is_hybrid=False)

    f0.TW.func.body = scalar
    f0.TW.___DO_push_all_to_instant___(0)
    f0.discretize()
    assert f0.error.L() < 0.01
    f1.TW.func.body = vector
    f1.TW.___DO_push_all_to_instant___(0)
    f1.discretize()
    assert f1.error.L() < 0.05
    f2.TW.func.body = vector
    f2.TW.___DO_push_all_to_instant___(0)
    f2.discretize()
    assert f2.error.L() < 0.07
    f3.TW.func.body = scalar
    f3.TW.___DO_push_all_to_instant___(0)
    f3.discretize()
    assert f3.error.L() < 0.05


    mesh = MeshGenerator('crazy', c=0.1)([5,5,5], EDM=None)
    space = SpaceInvoker('polynomials')([('Lobatto', 4), ('Lobatto', 4), ('Lobatto', 4)])
    FC = FormCaller(mesh, space)

    scalar = FC('scalar', p)
    vector = FC('vector', (u,v,w))
    f0 = FC('0-f', is_hybrid=False)
    f1 = FC('1-f', is_hybrid=False)
    f2 = FC('2-f', is_hybrid=False)
    f3 = FC('3-f', is_hybrid=False)

    f0.TW.func.body = scalar
    f0.TW.___DO_push_all_to_instant___(0)
    f0.discretize()
    assert f0.error.L() < 0.00005
    f1.TW.func.body = vector
    f1.TW.___DO_push_all_to_instant___(0)
    f1.discretize()
    assert f1.error.L() < 0.0008
    f2.TW.func.body = vector
    f2.TW.___DO_push_all_to_instant___(0)
    f2.discretize()
    assert f2.error.L() < 0.0008
    f3.TW.func.body = scalar
    f3.TW.___DO_push_all_to_instant___(0)
    f3.discretize()
    assert f3.error.L() < 0.0006

    return 1


def test_Form_NO1a_discretization_and_reconstruction():
    """"""
    if rAnk == mAster_rank:
        print("*** [test_Form_NO1a_discretization_and_reconstruction] ...... ", flush=True)

    def u(t, x, y, z): return np.cos(np.pi*x) + np.sin(np.pi*y) * np.sin(np.pi*z-0.125)**2 + t/2
    def v(t, x, y, z): return np.sin(np.pi*x) + np.sin(np.pi*y) * np.sin(np.pi*z-0.125)**2 + t/2
    def w(t, x, y, z): return np.sin(np.pi*x) + np.cos(np.pi*y) * np.cos(np.pi*z-0.125)**2 + t
    def p(t, x, y, z): return np.cos(np.pi*x) + np.sin(np.pi*y) * np.sin(np.pi*z-0.125)**2 + t/2

    try:
        mesh = MeshGenerator('LDC', l=1, w=1.2, h=1.5)([f'Lobatto:{3}', f'Lobatto:{4}', f'Lobatto:{5}'])
        space = SpaceInvoker('polynomials')([('Lobatto', 5), ('Lobatto', 5), ('Lobatto', 5)])
        FC = FormCaller(mesh, space)
        scalar = FC('scalar', p)
        vector = FC('vector', (u,v,w))
        f0 = FC('0-f', is_hybrid=False)
        f1 = FC('1-f', is_hybrid=False)
        f2 = FC('2-f', is_hybrid=False)
        f3 = FC('3-f', is_hybrid=False)
        f0.TW.func.body = scalar
        f0.TW.do.push_all_to_instant(0)
        f0.discretize()
        assert f0.error.L() < 0.00003
        f1.TW.func.body = vector
        f1.TW.do.push_all_to_instant(0)
        f1.discretize()
        assert f1.error.L() < 0.0004
        f2.TW.func.body = vector
        f2.TW.do.push_all_to_instant(0)
        f2.discretize()
        assert f2.error.L() < 0.0005
        f3.TW.func.body = scalar
        f3.TW.do.push_all_to_instant(0)
        f3.discretize()
        assert f3.error.L() < 0.0004
    except ThreeDimensionalTransfiniteInterpolationError:
        if rAnk == mAster_rank:
            print("    skip LDC mesh ...... ", flush=True)
        pass


    mesh = MeshGenerator('crazy_periodic', c=0.1)([3, 4, 5], EDM='chaotic')
    space = SpaceInvoker('polynomials')([('Lobatto', 5), ('Lobatto', 4), ('Lobatto', 3)])
    FC = FormCaller(mesh, space)

    scalar = FC('scalar', p)
    vector = FC('vector', (u,v,w))

    f0 = FC('0-f', is_hybrid=False)
    f1 = FC('1-f', is_hybrid=False)
    f2 = FC('2-f', is_hybrid=False)
    f3 = FC('3-f', is_hybrid=False)

    f0.TW.func.body = scalar
    f0.TW.do.push_all_to_instant(0)
    f0.discretize()
    assert f0.error.L() < 0.0004
    f1.TW.func.body = vector
    f1.TW.do.push_all_to_instant(0)
    f1.discretize()
    assert f1.error.L() < 0.004
    f2.TW.func.body = vector
    f2.TW.do.push_all_to_instant(0)
    f2.discretize()
    assert f2.error.L() < 0.0041
    f3.TW.func.body = scalar
    f3.TW.do.push_all_to_instant(0)
    f3.discretize()
    assert f3.error.L() < 0.003

    f0.cochain.globe = f0.cochain.globe
    f1.cochain.globe = f1.cochain.globe
    f2.cochain.globe = f2.cochain.globe
    f3.cochain.globe = f3.cochain.globe
    f0.TW.func.body = scalar
    f0.TW.do.push_all_to_instant(0)
    f0.discretize()
    assert f0.error.L() < 0.0004
    f1.TW.func.body = vector
    f1.TW.do.push_all_to_instant(0)
    f1.discretize()
    assert f1.error.L() < 0.004
    f2.TW.func.body = vector
    f2.TW.do.push_all_to_instant(0)
    f2.discretize()
    assert f2.error.L() < 0.0041
    f3.TW.func.body = scalar
    f3.TW.do.push_all_to_instant(0)
    f3.discretize()
    assert f3.error.L() < 0.003

    f0 = FC('0-f', is_hybrid=True)
    f1 = FC('1-f', is_hybrid=True)
    f2 = FC('2-f', is_hybrid=True)
    f3 = FC('3-f', is_hybrid=True)

    f0.TW.func.body = scalar
    f0.TW.do.push_all_to_instant(0)
    f0.discretize()
    assert f0.error.L() < 0.0004
    f1.TW.func.body = vector
    f1.TW.do.push_all_to_instant(0)
    f1.discretize()
    assert f1.error.L() < 0.004
    f2.TW.func.body = vector
    f2.TW.do.push_all_to_instant(0)
    f2.discretize()
    assert f2.error.L() < 0.0041
    f3.TW.func.body = scalar
    f3.TW.do.push_all_to_instant(0)
    f3.discretize()
    assert f3.error.L() < 0.003

    f0.cochain.globe = f0.cochain.globe
    f1.cochain.globe = f1.cochain.globe
    f2.cochain.globe = f2.cochain.globe
    f3.cochain.globe = f3.cochain.globe
    f0.TW.func.body = scalar
    f0.TW.do.push_all_to_instant(0)
    f0.discretize()
    assert f0.error.L() < 0.0004
    f1.TW.func.body = vector
    f1.TW.do.push_all_to_instant(0)
    f1.discretize()
    assert f1.error.L() < 0.004
    f2.TW.func.body = vector
    f2.TW.do.push_all_to_instant(0)
    f2.discretize()
    assert f2.error.L() < 0.0041
    f3.TW.func.body = scalar
    f3.TW.do.push_all_to_instant(0)
    f3.discretize()
    assert f3.error.L() < 0.003


    mesh = MeshGenerator('bridge_arch_cracked',)([3,2,4], EDM='chaotic')
    space = SpaceInvoker('polynomials')([('Lobatto', 3), ('Lobatto', 4), ('Lobatto', 2)])
    FC = FormCaller(mesh, space)

    scalar = FC('scalar', p)
    vector = FC('vector', (u,v,w))
    f0 = FC('0-f', is_hybrid=False)
    f1 = FC('1-f', is_hybrid=False)
    f2 = FC('2-f', is_hybrid=False)
    f3 = FC('3-f', is_hybrid=False)

    f0.TW.func.body = scalar
    f0.TW.do.push_all_to_instant(0)
    f0.discretize()
    assert f0.error.L() < 0.01
    f1.TW.func.body = vector
    f1.TW.do.push_all_to_instant(0)
    f1.discretize()
    assert f1.error.L() < 0.05
    f2.TW.func.body = vector
    f2.TW.do.push_all_to_instant(0)
    f2.discretize()
    assert f2.error.L() < 0.07
    f3.TW.func.body = scalar
    f3.TW.do.push_all_to_instant(0)
    f3.discretize()
    assert f3.error.L() < 0.05

    return 1


def test_Form_NO1b_trace_form_Rd_and_Rc():
    """"""
    if rAnk == mAster_rank:
        print("*** [test_Form_NO1b_trace_form_Rd_and_Rc] ...... ", flush=True)

    def p(t, x, y, z):
        return np.cos(2*np.pi*x) + np.cos(np.pi*y) * np.cos(np.pi*z-0.125)**2 + t/2

    mesh = MeshGenerator('crazy', c=0.0)([3, 4, 5])
    space = SpaceInvoker('polynomials')([('Lobatto', 9),
                                         ('Lobatto', 8),
                                         ('Lobatto', 7)])
    FC = FormCaller(mesh, space)
    flux = FC('scalar', p)
    t = 0
    # test 0-Trace form reconstruct ...............
    t0 = FC('0-t')
    t0.TW.func.do.set_func_body_as(flux)
    t0.TW.current_time = t
    t0.TW.do.push_all_to_instant()
    t0.discretize()
    xi = eta = sigma = np.linspace(-1, 1, 50)
    xyz, V = t0.reconstruct(xi, eta, sigma)
    for i in xyz:
        x, y, z = xyz[i]
        v = V[i][0]
        v_exact = p(t, x, y, z)
        assert np.max(np.abs(v - v_exact)) < 1e-5 # must be accurate enough!


    # test 2-Trace form reconstruct...............
    t2 = FC('2-t')
    t2.TW.func.do.set_func_body_as(flux)
    t2.TW.current_time = t
    t2.TW.do.push_all_to_instant()
    t2.discretize()
    xi = eta = sigma = np.linspace(-1, 1, 50)
    xyz, V = t2.reconstruct(xi, eta, sigma)
    for i in xyz:
        x, y, z = xyz[i]
        v = V[i][0]
        v_exact = p(t, x, y, z)
        assert np.max(np.abs(v - v_exact)) < 1e-4 # must be accurate enough!

    # ------------------------------------------------------------------------------------------

    def uuu(t, x, y, z):
        return 5*np.sin(0.89*np.pi*x) + 6*np.cos(np.pi*y) + 7*np.cos(np.pi*z-0.578)**2 + t*8.5
    def vvv(t, x, y, z):
        return 5*np.cos(2.21*np.pi*x) + 6*np.sin(np.pi*y) + 7*np.cos(np.pi*z-0.12)**2 + t*8
    def www(t, x, y, z):
        return 5*np.cos(np.pi*x) + 6*np.cos(np.pi*y) + 7*np.sin(np.pi*z-0.15)**2 + t*10

    if rAnk == mAster_rank:
        load = random.randint(500,1000)
        t = random.random()
    else:
        load, t = None, None
    load, t = cOmm.bcast([load, t], root=mAster_rank)
    mesh, space = random_3D_mesh_and_space_of_total_load_around(load, exclude_periodic = True)
    FC = FormCaller(mesh, space)
    flux = FC('scalar', p)
    velo = FC('vector', (uuu, vvv, www))

    t0 = FC('0-t')
    t1 = FC('1-t')
    t2 = FC('2-t')

    f0 = FC('0-f', is_hybrid=True)
    f1 = FC('1-f', is_hybrid=True)
    f2 = FC('2-f', is_hybrid=True)

    S0 = t0.matrices.selective
    S1 = t1.matrices.selective
    S2 = t2.matrices.selective

    # t0 & f0, discretization and selective matrix
    t0.TW.func.do.set_func_body_as(flux)
    t0.TW.current_time = t
    t0.TW.do.push_all_to_instant()
    t0.discretize()
    f0.TW.func.do.set_func_body_as(flux)
    f0.TW.current_time = t
    f0.TW.do.push_all_to_instant()
    f0.discretize()
    t0_local = t0.cochain.local
    f0_local = f0.cochain.local
    for i in t0_local:
        A0 = t0_local[i] -  S0[i] @ f0_local[i]
        np.testing.assert_array_almost_equal(A0, 0)

    # t1 & f1, discretization and selective matrix
    t1.TW.func.do.set_func_body_as(velo)
    t1.TW.current_time = t
    t1.TW.do.push_all_to_instant()
    t1.discretize()
    f1.TW.func.do.set_func_body_as(velo)
    f1.TW.current_time = t
    f1.TW.do.push_all_to_instant()
    f1.discretize()
    t1_local = t1.cochain.local
    f1_local = f1.cochain.local
    for i in t1_local:
        A1 = t1_local[i] -  S1[i] @ f1_local[i]
        np.testing.assert_array_almost_equal(A1, 0)

    # t2 & f2, discretization and selective matrix
    t2.TW.func.do.set_func_body_as(velo)
    t2.TW.current_time = t
    t2.TW.do.push_all_to_instant()
    t2.discretize()
    f2.TW.func.do.set_func_body_as(velo)
    f2.TW.current_time = t
    f2.TW.do.push_all_to_instant()
    f2.discretize()
    t2_local = t2.cochain.local
    f2_local = f2.cochain.local
    for i in t2_local:
        A2 = t2_local[i] -  S2[i] @ f2_local[i]
        np.testing.assert_array_almost_equal(A2, 0)

    return 1


def test_Form_NO2_mass_matrix():
    """
    Unittests for the mesh.
    """
    if rAnk == mAster_rank:
        print("*** [test_Form_NO2_mass_matrix] ...... ", flush=True)

    mesh = MeshGenerator('crazy', c=0.24)([2, 2, 2], EDM='debug')
    space = SpaceInvoker('polynomials')([('Lobatto', 2), ('Lobatto', 2), ('Lobatto', 2)])
    FC = FormCaller(mesh, space)

    f0 = FC('0-f', is_hybrid=False)
    f1 = FC('1-f', is_hybrid=False)
    f2 = FC('2-f', is_hybrid=True)
    f3 = FC('3-f', is_hybrid=True)

    benchmark0 = np.array([
        3.69464136e-04,  1.88796948e-04, -8.62687140e-05,  1.88796948e-04,
        9.64309140e-05, -4.41505770e-05, -8.62687140e-05, -4.41505770e-05,
        2.00428485e-05,  1.88796948e-04,  9.64309140e-05, -4.41505770e-05,
        9.64309140e-05,  4.92316770e-05, -2.25833985e-05, -4.41505770e-05,
       -2.25833985e-05,  1.02754793e-05, -8.62687140e-05, -4.41505770e-05,
        2.00428485e-05, -4.41505770e-05, -2.25833985e-05,  1.02754793e-05,
        2.00428485e-05,  1.02754793e-05, -4.62962963e-06])
    benchmark1 = np.array([
        2.46972792e-02, -2.03100085e-03,  1.27185565e-02, -9.65220878e-04,
       -5.61944452e-03,  5.83169531e-04,  1.27185565e-02, -9.65220878e-04,
        6.54630854e-03, -4.55572897e-04, -2.89909366e-03,  2.81861533e-04,
       -5.61944452e-03,  5.83169531e-04, -2.89909366e-03,  2.81861533e-04,
        1.27080404e-03, -1.60377232e-04, -2.39174481e-03, -3.66693902e-03,
        2.18959883e-04,  4.66797763e-04,  1.72673100e-04,  3.75541827e-05,
       -1.19958574e-03, -1.88465494e-03,  1.10994391e-04,  2.34794758e-04,
        1.13900041e-04,  1.81982257e-05,  5.92366200e-04,  8.39956606e-04,
       -5.24682962e-05, -1.14605626e-04, -1.82303885e-06, -1.02568442e-05,
       -2.39174481e-03, -3.66693902e-03,  2.18959883e-04, -1.19958574e-03,
       -1.88465494e-03,  1.10994391e-04,  5.92366200e-04,  8.39956606e-04,
       -5.24682962e-05,  4.66797763e-04,  1.72673100e-04,  3.75541827e-05,
        2.34794758e-04,  1.13900041e-04,  1.81982257e-05, -1.14605626e-04,
       -1.82303885e-06, -1.02568442e-05])
    benchmark2 = np.array([
        1.54538742,  0.79471961, -0.35330801, -0.1544278 , -0.08356559,
        0.02907942, -0.1544278 , -0.08356559,  0.02907942,  0.0447084 ,
        0.0129169 , -0.02533305,  0.22911758, -0.04166723,  0.34090232,
        0.0078891 , -0.01563826, -0.00465875,  0.0323816 , -0.00314733,
        0.02466658,  0.03731131,  0.00258758, -0.00262363,  0.22911758,
       -0.04166723,  0.0323816 , -0.00314733,  0.34090232,  0.0078891 ,
        0.02466658,  0.03731131, -0.01563826, -0.00465875,  0.00258758,
       -0.00262363])
    benchmark3 = np.array([
        7.76343952e+01, -1.70094465e+01, -1.70094465e+01, -2.48764637e-02,
       -1.70094465e+01, -2.48764637e-02, -2.48764637e-02, -3.05021211e+00])

    M0 = f0.matrices.mass
    M1 = f1.matrices.mass
    M2 = f2.matrices.mass
    M3 = f3.matrices.mass
    if 0 in mesh.elements:
        np.testing.assert_array_almost_equal(M0[0].toarray()[0,:], benchmark0)
        np.testing.assert_array_almost_equal(M1[0].toarray()[0,:], benchmark1)
        np.testing.assert_array_almost_equal(M2[0].toarray()[0,:], benchmark2)
        np.testing.assert_array_almost_equal(M3[0].toarray()[0,:], benchmark3)

    M0_3 = M0 * 3
    three_M1 = 3 * M1
    MMM = M2*5 - 2*M2
    M3M5 = M3*5 + 2*M3
    M13 = M1 / 3
    mM14 = - M1 / 4
    if 0 in mesh.elements:
        np.testing.assert_array_almost_equal(M0_3[0].toarray()[0,:], 3*benchmark0)
        np.testing.assert_array_almost_equal(three_M1[0].toarray()[0,:], 3*benchmark1)
        np.testing.assert_array_almost_equal(MMM[0].toarray()[0,:], 3*benchmark2)
        np.testing.assert_array_almost_equal(M3M5[0].toarray()[0,:], 7*benchmark3)
        np.testing.assert_array_almost_equal(M13[0].toarray()[0,:], benchmark1/3)
        np.testing.assert_array_almost_equal(mM14[0].toarray()[0,:], -0.25*benchmark1)

    M2T = M2.T
    assert M2.gathering_matrices == (M2T.gathering_matrices[1], M2T.gathering_matrices[0])
    assert M0_3.gathering_matrices == M0.gathering_matrices
    assert three_M1.gathering_matrices == M1.gathering_matrices == M13.gathering_matrices == mM14.gathering_matrices

    if 0 in mesh.elements:
        np.testing.assert_array_almost_equal(M2[0].toarray()[:,0], benchmark2)

    mM2T = - M2.T
    assert mM2T.gathering_matrices == M2T.gathering_matrices
    if 0 in mesh.elements:
        np.testing.assert_array_almost_equal(mM2T[0].toarray()[:,0], -benchmark2)

    return 1


def test_Form_NO3_incidence_matrices():
    """
    Unittests for the mesh.
    """
    if rAnk == mAster_rank:
        print("*** [test_Form_NO3_incidence_matrices] ...... ", flush=True)

    mesh = MeshGenerator('crazy', c=0.0)([2, 3, 4], EDM='debug')
    space = SpaceInvoker('polynomials')([('Lobatto', 8), ('Lobatto', 7), ('Lobatto', 6)])
    FC = FormCaller(mesh, space)
    es = ExactSolutionSelector(mesh)('icpsNS:sincosRD')
    np.testing.assert_almost_equal(es.status.helicity(0), -12.566370614359162)
    np.testing.assert_almost_equal(es.status.kinetic_energy(0), 1.4999999999999998)
    np.testing.assert_almost_equal(es.status.enstrophy(0), 59.217626406536155)
    # so kinetic_energy = 0.5 * (L^2 norm of velocity)**2
    np.testing.assert_almost_equal(
        es.status.___PRIVATE_compute_L_norm_of___('velocity', time=0, n=2), np.sqrt(3))

    f0 = FC('0-f', is_hybrid=True)

    f0.TW.func.body = es.status.pressure
    f0.TW.___DO_push_all_to_instant___(0)
    f0.discretize()
    assert f0.error.L() < 5.8e-7

    f1 = f0.coboundary()
    f1.TW.func.body = es.status.gradient_of_pressure
    f1.TW.___DO_push_all_to_instant___(0)
    assert f1.error.L() < 3.6e-5

    f1.TW.func.body = es.status.velocity
    f1.TW.___DO_push_all_to_instant___(0)
    f1.discretize()
    f2 = f1.coboundary()
    f2.TW.func.body = es.status.vorticity
    f2.TW.___DO_push_all_to_instant___(0)
    assert f2.error.L() < 6.2e-5

    f3 = f2.coboundary()
    for i in mesh.elements:
        assert np.max(np.abs(f3.cochain.local[i])) < 1e-15

    E10 = f0.coboundary.incidence_matrix
    E21 = f1.coboundary.incidence_matrix
    E32 = f2.coboundary.incidence_matrix

    E21E10 = E21 @ E10
    E32E21 = E32 @ E21
    for i in mesh.elements:
        assert E21E10[i].nnz == 0
        assert E32E21[i].nnz == 0

    f3.TW.func.body = es.status.divergence_of_velocity
    f3.TW.do.push_all_to_instant(0)
    L_inf = f3.error.L(n='infinity')
    assert L_inf < 1e-8

    # we test curl of vorticity for incompressible NS ----------- BELOW ------------------------------------------------
    f1 = FC('1-f', is_hybrid=False)
    f1.TW.func.body = es.status.vorticity
    f1.TW.do.push_all_to_instant(0)
    f1.discretize()
    f2 = f1.coboundary()
    f2.TW.func.body = es.status.curl_of_vorticity
    f2.TW.do.push_all_to_instant(0)
    assert f2.error.L() < 0.00045

    mesh = MeshGenerator('crazy', c=0.0)([6, 6, 6], EDM='debug')
    space = SpaceInvoker('polynomials')([('Lobatto', 8), ('Lobatto', 8), ('Lobatto', 8)])
    FC = FormCaller(mesh, space)
    es = ExactSolutionSelector(mesh)('icpsNS:sincosRD')

    f1 = FC('1-f', is_hybrid=True)
    f1.TW.func.body = es.status.vorticity
    f1.TW.do.push_all_to_instant(0)
    f1.discretize()
    f2 = f1.coboundary()
    f2.TW.func.body = es.status.curl_of_vorticity
    f2.TW.do.push_all_to_instant(0)
    assert f2.error.L() < 1e-7

    # we test numerical gradient, curl and div with coboundary of standard forms------------------------------------
    if rAnk == mAster_rank:
        c = random.random() / 10
        if c < 0.05:
            c = 0
        else:
            pass
    else:
        c = None
    c = cOmm.bcast(c, root=mAster_rank)
    t = 0

    mesh = MeshGenerator('crazy', c=c)([5, 4, 5], EDM='debug')
    space = SpaceInvoker('polynomials')([('Lobatto', 4), ('Lobatto', 5), ('Lobatto', 4)])
    FC = FormCaller(mesh, space)

    def w0(t, x, y, z): return -2*np.pi * np.sin(x) * np.cos(y) * np.cos(z*2.45) + np.sin(t)
    def w1(t, x, y, z): return -3*np.pi * np.cos(x/np.pi) * np.sin(0.5*y) * np.cos(z/2) + np.sin(2*t)
    def w2(t, x, y, z): return -np.pi * np.cos(1.11*x) * np.sin(y) * np.sin(np.pi*z) + (1.5-np.cos(t/2))
    def p(t, x, y, z): return 3**np.pi * np.sin(2*x) * np.sin(y/3) * np.sin(z*3.45) + np.sin(0.57*t)

    V = FC('vector', (w0, w1, w2))
    S = FC('scalar', p)
    gS = S.numerical.gradient
    cV = V.numerical.curl
    dV = V.numerical.divergence

    f0 = FC('0-f', is_hybrid=False)
    f0.TW.func.body = S
    f0.TW.do.push_all_to_instant(t)
    f0.discretize()
    f1 = f0.coboundary()
    f1.TW.func.body = gS
    f1.TW.do.push_all_to_instant(t)
    assert f1.error.L() < 0.0065

    f1 = FC('1-f', is_hybrid=False)
    f1.TW.func.body = V
    f1.TW.do.push_all_to_instant(t)
    f1.discretize()
    f2 = f1.coboundary()
    f2.TW.func.body = cV
    f2.TW.do.push_all_to_instant(t)
    assert f2.error.L() < 0.002

    f2 = FC('2-f', is_hybrid=False)
    f2.TW.func.body = V
    f2.TW.do.push_all_to_instant(t)
    f2.discretize()
    f3 = f2.coboundary()
    f3.TW.func.body = dV
    f3.TW.do.push_all_to_instant(t)
    assert f3.error.L() < 0.0025

    # +++++++++++++++++++++++++++++++++++++++++ ABOVE ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    return 1


def test_Form_NO4_cross_product_1():
    """
    Unittests for the special method ``cross_product`` of the standard 1-form.
    """
    if rAnk == mAster_rank:
        print("*** [test_Form_NO4_cross_product_1] ...... ", flush=True)

    mesh = MeshGenerator('crazy', c=0.1)([3, 3, 4], EDM='debug')
    space = SpaceInvoker('polynomials')([('Lobatto', 3), ('Lobatto', 3), ('Lobatto', 2)])
    fmCa = FormCaller(mesh, space)

    def w0(t, x, y, z): return -np.pi * np.sin(x) * np.cos(y) * np.cos(z) + t
    def w1(t, x, y, z): return -np.pi * np.cos(x) * np.sin(y) * np.cos(z) + t
    def w2(t, x, y, z): return -np.pi * np.cos(x) * np.sin(y) * np.sin(z) + t
    def u0(t, x, y, z): return np.sin(np.pi*x) * np.sin(np.pi*y) * np.cos(np.pi*z) + t
    def u1(t, x, y, z): return np.cos(np.pi*x) * np.sin(np.pi*y) * np.sin(np.pi*z) + t
    def u2(t, x, y, z): return np.sin(np.pi*x) * np.cos(np.pi*y) * np.sin(np.pi*z) + t
    # x = w X u
    def x0(t, x, y, z): return w1(t, x, y, z) * u2(t, x, y, z) - w2(t, x, y, z) * u1(t, x, y, z)
    def x1(t, x, y, z): return w2(t, x, y, z) * u0(t, x, y, z) - w0(t, x, y, z) * u2(t, x, y, z)
    def x2(t, x, y, z): return w0(t, x, y, z) * u1(t, x, y, z) - w1(t, x, y, z) * u0(t, x, y, z)

    vw = fmCa('vector', func=(w0, w1, w2))
    vu = fmCa('vector', func=(u0, u1, u2))
    vx = fmCa('vector', func=(x0, x1, x2))

    u = fmCa('1-f')
    w = fmCa('1-f')
    # x:= w times u
    x = fmCa('1-f')

    w.TW.func.body = vw
    w.TW.___DO_push_all_to_instant___(0)
    w.discretize()

    # represent (w times u, e) or (w^1 wedge u^1, e^1)
    MW = w.special.cross_product_1f__ip_1f(u, x)
    M = x.matrices.mass

    u.TW.func.body = vu
    u.TW.___DO_push_all_to_instant___(0)
    u.discretize()
    _xcl = {}
    for i in mesh.elements:
        P = spspalinalg.inv(M[i].tocsc()) @ MW[i]
        _xcl[i] = P @ u.cochain.local[i]
    x.cochain.local = _xcl

    x.TW.func.body = vx
    x.TW.___DO_push_all_to_instant___(0)
    assert x.error.L() < 0.04

    # now, test the orthogonal mesh ...
    mesh = MeshGenerator('crazy', c=0.0)([3, 3, 4], EDM='debug')
    fmCa = FormCaller(mesh, space)
    u = fmCa('1-f')
    w = fmCa('1-f')
    # # x:= w times u
    x = fmCa('1-f')

    vw = fmCa('vector', func=(w0, w1, w2))
    vu = fmCa('vector', func=(u0, u1, u2))
    vx = fmCa('vector', func=(x0, x1, x2))

    w.TW.func.body = vw
    w.TW.___DO_push_all_to_instant___(0)
    w.discretize()

    # represent (w times u, e) or (w^1 wedge u^1, e^1)
    MW = w.special.cross_product_1f__ip_1f(u, x)
    M = x.matrices.mass

    u.TW.func.body = vu
    u.TW.___DO_push_all_to_instant___(0)
    u.discretize()
    _xcl = {}
    if mesh.elements.num > 0:
        i = mesh.elements.indices[0]
        invM = spspalinalg.inv(M[i].tocsc())

    for i in mesh.elements:
        # noinspection PyUnboundLocalVariable
        P = invM @ MW[i]
        _xcl[i] = P @ u.cochain.local[i]

    x.cochain.local = _xcl

    x.TW.func.body = vx
    x.TW.___DO_push_all_to_instant___(0)
    assert x.error.L() < 0.022

    return 1


def test_Form_NO5_cross_product_2():
    """
    Unittests for the special method ``cross_product`` of the standard 1-form.
    """
    if rAnk == mAster_rank:
        print("*** [test_Form_NO5_cross_product_2] ...... ", flush=True)

    mesh = MeshGenerator('crazy', c=0.0, bounds=([0,1],[0,1],[0,1]))([2,3,1], EDM='debug')
    space = SpaceInvoker('polynomials')([('Lobatto', 5), ('Lobatto', 4), ('Lobatto', 6)])
    fmCa = FormCaller(mesh, space)

    def w0(t, x, y, z): return -np.pi * np.sin(x) * np.cos(y) * np.cos(z) + t
    def w1(t, x, y, z): return -np.pi * np.cos(x) * np.sin(y) * np.cos(z) + t
    def w2(t, x, y, z): return -np.pi * np.cos(x) * np.sin(y) * np.sin(z) + t
    def u0(t, x, y, z): return np.sin(np.pi*x) * np.sin(y) * np.cos(np.pi*z) + t
    def u1(t, x, y, z): return np.cos(x) * np.sin(np.pi*y) * np.sin(np.pi*z) + t
    def u2(t, x, y, z): return np.sin(np.pi*x) * np.cos(np.pi*y) * np.sin(z) + t
    # x = w X u
    def x0(t, x, y, z): return w1(t, x, y, z) * u2(t, x, y, z) - w2(t, x, y, z) * u1(t, x, y, z)
    def x1(t, x, y, z): return w2(t, x, y, z) * u0(t, x, y, z) - w0(t, x, y, z) * u2(t, x, y, z)
    def x2(t, x, y, z): return w0(t, x, y, z) * u1(t, x, y, z) - w1(t, x, y, z) * u0(t, x, y, z)

    vw = fmCa('vector', func=(w0, w1, w2))
    vu = fmCa('vector', func=(u0, u1, u2))
    vx = fmCa('vector', func=(x0, x1, x2))

    w = fmCa('2-f')
    u = fmCa('2-f')
    # x:= w times u
    x = fmCa('2-f')

    w.TW.func.body = vw
    w.TW.___DO_push_all_to_instant___(0)
    w.discretize()

    # represent (w times u, e) or (star w^2 wedge star u^2, e^2)
    MW = w.special.cross_product_2f__ip_2f(u, x)
    M = x.matrices.mass

    u.TW.func.body = vu
    u.TW.___DO_push_all_to_instant___(0)
    u.discretize()

    _xcl = {}

    if mesh.elements.num > 0:
        i = mesh.elements.indices[0]
        invM = spspalinalg.inv(M[i].tocsc())

    for i in mesh.elements:
        # noinspection PyUnboundLocalVariable
        P = invM @ MW[i]
        _xcl[i] = P @ u.cochain.local[i]
    x.cochain.local = _xcl

    x.TW.func.body = vx
    x.TW.___DO_push_all_to_instant___(0)
    assert x.error.L() < 0.00085

    return 1


def test_Form_NOx_cross_product_3():
    """
    Unittests for the special method ``cross_product`` of the standard 1-form.
    """
    if rAnk == mAster_rank:
        print("*** [test_Form_NOx_cross_product_3] ...... ", flush=True)

    mesh = MeshGenerator('crazy', c=0.0, bounds=([0,1],[0,1],[0,1]))([10,10,10], EDM=None)
    space = SpaceInvoker('polynomials')([('Lobatto', 2), ('Lobatto', 2), ('Lobatto', 2)])
    FC = FormCaller(mesh, space)

    def w0(t, x, y, z): return -np.pi * np.sin(x) * np.cos(y) * np.cos(z) + t
    def w1(t, x, y, z): return -np.pi * np.cos(x) * np.sin(y) * np.cos(z) + t
    def w2(t, x, y, z): return -np.pi * np.cos(x) * np.sin(y) * np.sin(z) + t
    def u0(t, x, y, z): return np.sin(np.pi*x) * np.sin(y) * np.cos(np.pi*z) + t
    def u1(t, x, y, z): return np.cos(x) * np.sin(np.pi*y) * np.sin(np.pi*z) + t
    def u2(t, x, y, z): return np.sin(np.pi*x) * np.cos(np.pi*y) * np.sin(z) + t
    # x = w X u
    def x0(t, x, y, z): return w1(t, x, y, z) * u2(t, x, y, z) - w2(t, x, y, z) * u1(t, x, y, z)
    def x1(t, x, y, z): return w2(t, x, y, z) * u0(t, x, y, z) - w0(t, x, y, z) * u2(t, x, y, z)
    def x2(t, x, y, z): return w0(t, x, y, z) * u1(t, x, y, z) - w1(t, x, y, z) * u0(t, x, y, z)

    vw = FC('vector', func=(w0, w1, w2))
    vu = FC('vector', func=(u0, u1, u2))
    vx = FC('vector', func=(x0, x1, x2))

    w = FC('1-f')
    u = FC('2-f')
    # x:= w times u
    x = FC('2-f')

    w.TW.func.body = vw
    w.TW.do.push_all_to_instant(0)
    w.discretize()

    # represent (w1 times u2, e2)
    MW = w.special.cross_product_2f__ip_2f(u, x)
    M = x.matrices.mass

    u.TW.func.body = vu
    u.TW.do.push_all_to_instant(0)
    u.discretize()

    _xcl = dict()

    if mesh.elements.num > 0:
        i = mesh.elements.indices[0]
        invM = spspalinalg.inv(M[i].tocsc())

    for i in mesh.elements:
        # noinspection PyUnboundLocalVariable
        P = invM @ MW[i]
        _xcl[i] = P @ u.cochain.local[i]

        result = np.einsum('j, jk, k ->', u.cochain.local[i], MW[i].toarray(), u.cochain.local[i],
                           optimize='greedy')
        np.testing.assert_almost_equal(result, 0)


    x.cochain.local = _xcl

    x.TW.func.body = vx
    x.TW.do.push_all_to_instant(0)
    assert x.error.L() < 0.0055

    return 1


def test_Form_NO6_resemble():
    """
    Unittests for the special method ``DO_imitate`` of the standard 1-form.
    """
    if rAnk == mAster_rank:
        print("*** [test_Form_NO6_resemble] ...... ", flush=True)
    def p(t, x, y, z):
        return np.cos(2*np.pi*x) * np.cos(2*np.pi*y) * np.cos(2*np.pi*z) + t/2
    def u(t, x, y, z):
        return np.sin(2*np.pi*x) * np.cos(2*np.pi*y) * np.cos(2*np.pi*z) + t/2
    def v(t, x, y, z):
        return np.cos(2*np.pi*x) * np.sin(2*np.pi*y) * np.cos(2*np.pi*z) + t/2
    def w(t, x, y, z):
        return np.cos(2*np.pi*x) * np.cos(2*np.pi*y) * np.sin(2*np.pi*z) + t/2

    mesh = MeshGenerator('crazy', c=0.0)([2,2,2], EDM='debug')
    space = SpaceInvoker('polynomials')([('Lobatto',4), ('Lobatto',4), ('Lobatto',4)])
    fc = FormCaller(mesh, space)

    scalar = fc('scalar', p)
    vector = fc('vector', (u,v,w))
    f0 = fc('0-f', is_hybrid=False)
    f1 = fc('1-f', is_hybrid=False)
    f2 = fc('2-f', is_hybrid=True)
    f3 = fc('3-f', is_hybrid=True)

    f0.TW.func.body = scalar
    f0.TW.do.push_all_to_instant(0)
    f0.discretize()
    f1.TW.func.body = vector
    f1.TW.do.push_all_to_instant(0)
    f1.discretize()
    f2.TW.func.body = vector
    f2.TW.do.push_all_to_instant(0)
    f2.discretize()
    f3.TW.func.body = scalar
    f3.TW.do.push_all_to_instant(0)
    f3.discretize()

    save(f0, 'f0')
    save(f1, 'f1')
    save(f2, 'f2')
    save(f3, 'f3')
    l0 = f0.do.compute_L2_energy_with()
    l1 = f1.do.compute_L2_energy_with()
    l2 = f2.do.compute_L2_energy_with()
    l3 = f3.do.compute_L2_energy_with()

    F0 = fc('0-f', is_hybrid=True)
    F1 = fc('1-f', is_hybrid=True)
    F2 = fc('2-f', is_hybrid=False)
    F3 = fc('3-f', is_hybrid=False)
    F0.do.resemble(f0)
    F1.do.resemble(f1)
    F2.do.resemble(f2)
    F3.do.resemble(f3)
    np.testing.assert_almost_equal(F0.do.compute_Ln_diff_from(f0), 0)
    np.testing.assert_almost_equal(F1.do.compute_Ln_diff_from(f1), 0)
    np.testing.assert_almost_equal(F2.do.compute_Ln_diff_from(f2), 0)
    np.testing.assert_almost_equal(F3.do.compute_Ln_diff_from(f3), 0)

    mesh = MeshGenerator('crazy', c=0.0)([3,2,3], EDM='debug')
    space = SpaceInvoker('polynomials')([('Lobatto',4), ('Lobatto',5), ('Lobatto',4)])
    FC = FormCaller(mesh, space)
    F0 = FC('0-f', is_hybrid=False)
    F1 = FC('1-f', is_hybrid=False)
    F2 = FC('2-f', is_hybrid=True)
    F3 = FC('3-f', is_hybrid=False)
    F0.do.resemble(f0)
    F1.do.resemble(f1)
    F2.do.resemble(f2, density=20000)
    F3.do.resemble(f3, density=20000)
    L0 = F0.do.compute_L2_energy_with()
    L1 = F1.do.compute_L2_energy_with()
    L2 = F2.do.compute_L2_energy_with()
    L3 = F3.do.compute_L2_energy_with()
    assert (np.abs(l0-L0)/l0) < 0.03
    assert (np.abs(l1-L1)/l1) < 0.02
    assert (np.abs(l2-L2)/l2) < 0.01
    assert (np.abs(l3-L3)/l3) < 0.03
    F0 = FC('0-f', is_hybrid=False)
    F1 = FC('1-f', is_hybrid=False)
    F2 = FC('2-f', is_hybrid=True)
    F3 = FC('3-f', is_hybrid=False)
    F0.do.resemble('f0.mi')
    F1.do.resemble('f1.mi')
    F2.do.resemble('f2.mi', density=20000)
    F3.do.resemble('f3.mi', density=20000)
    L0 = F0.do.compute_L2_energy_with()
    L1 = F1.do.compute_L2_energy_with()
    L2 = F2.do.compute_L2_energy_with()
    L3 = F3.do.compute_L2_energy_with()
    assert (np.abs(l0-L0)/l0) < 0.03
    assert (np.abs(l1-L1)/l1) < 0.02
    assert (np.abs(l2-L2)/l2) < 0.01
    assert (np.abs(l3-L3)/l3) < 0.03

    t2 = FC('2-t')
    t2.TW.func.body = scalar
    t2.TW.do.push_all_to_instant(0)
    t2.discretize()
    save(t2, 'trace')
    T2 = FC('2-t')
    T2.do.resemble(t2)
    for i in mesh.elements:
        assert np.abs(np.sum(T2.cochain.local[i]- t2.cochain.local[i])) < 0.002
    T2 = FC('2-t')
    T2.do.resemble('trace.mi')
    for i in mesh.elements:
        assert np.abs(np.sum(T2.cochain.local[i]- t2.cochain.local[i])) < 0.002

    if rAnk == mAster_rank:
        os.remove('trace.mi')
        os.remove('f0.mi')
        os.remove('f1.mi')
        os.remove('f2.mi')
        os.remove('f3.mi')

    return 1


def test_Form_No7_with_other_element_numbering_AUTO():
    if rAnk == mAster_rank:
        print("*** [test_Form_No7_with_other_element_numbering---AUTO] ...... ", flush=True)

    def u(t, x, y, z): return np.cos(np.pi*x) + np.sin(np.pi*y) * np.sin(np.pi*z-0.125)**2 + t/2
    def v(t, x, y, z): return np.sin(np.pi*x) + np.sin(np.pi*y) * np.sin(np.pi*z-0.125)**2 + t/2
    def w(t, x, y, z): return np.sin(np.pi*x) + np.cos(np.pi*y) * np.cos(np.pi*z-0.125)**2 + t
    def p(t, x, y, z): return np.cos(np.pi*x) + np.sin(np.pi*y) * np.sin(np.pi*z-0.125)**2 + t/2

    if rAnk == mAster_rank:
        t = random.uniform(-100,100)
    else:
        t = None
    t = cOmm.bcast(t, root=mAster_rank)

    if rAnk == mAster_rank:
        i = random.randint(2,4)
        j = random.randint(2,4)
        k = random.randint(2,4)
    else:
        i, j, k = None, None, None
    i, j, k = cOmm.bcast([i, j, k], root=mAster_rank)
    space = SpaceInvoker('polynomials')([('Lobatto', i), ('Lobatto', j), ('Lobatto', k)])

    if rAnk == mAster_rank:
        i = random.randint(5,8)
        j = random.randint(5,9)
        k = random.randint(6,7)
    else:
        i, j, k = None, None, None
    i, j, k = cOmm.bcast([i, j, k], root=mAster_rank)
    MESH = MeshGenerator('bridge_arch_cracked')([i, j, k], EDM='debug', show_info=False)
    mesh = MeshGenerator('bridge_arch_cracked')([i, j, k], EDM=None, show_info=False)
    MQ = MESH.___PRIVATE_element_division_and_numbering_quality___()[0]
    mQ = mesh.___PRIVATE_element_division_and_numbering_quality___()[0]
    assert MQ <= mQ, "We expect this!"

    FC = FormCaller(mesh, space)
    scalar = FC('scalar', p)
    vector = FC('vector', (u,v,w))
    f0 = FC('0-f', is_hybrid=False)
    f1 = FC('1-f', is_hybrid=False)
    f2 = FC('2-f', is_hybrid=False)
    f3 = FC('3-f', is_hybrid=False)
    f0.TW.func.body = scalar
    f0.TW.___DO_push_all_to_instant___(t)
    f0.discretize()
    F0E = f0.error.L()
    f1.TW.func.body = vector
    f1.TW.___DO_push_all_to_instant___(t)
    f1.discretize()
    F1E = f1.error.L()
    f2.TW.func.body = vector
    f2.TW.___DO_push_all_to_instant___(t)
    f2.discretize()
    F2E = f2.error.L()
    f3.TW.func.body = scalar
    f3.TW.___DO_push_all_to_instant___(t)
    f3.discretize()
    F3E = f3.error.L()


    FC = FormCaller(MESH, space)
    scalar = FC('scalar', p)
    vector = FC('vector', (u,v,w))
    f0 = FC('0-f', is_hybrid=False)
    f1 = FC('1-f', is_hybrid=False)
    f2 = FC('2-f', is_hybrid=False)
    f3 = FC('3-f', is_hybrid=False)
    f0.TW.func.body = scalar
    f0.TW.___DO_push_all_to_instant___(t)
    f0.discretize()
    F0e = f0.error.L()
    f1.TW.func.body = vector
    f1.TW.___DO_push_all_to_instant___(t)
    f1.discretize()
    F1e = f1.error.L()
    f2.TW.func.body = vector
    f2.TW.___DO_push_all_to_instant___(t)
    f2.discretize()
    F2e = f2.error.L()
    f3.TW.func.body = scalar
    f3.TW.___DO_push_all_to_instant___(t)
    f3.discretize()
    F3e = f3.error.L()

    np.testing.assert_almost_equal(F0E, F0e)
    np.testing.assert_almost_equal(F1E, F1e)
    np.testing.assert_almost_equal(F2E, F2e)
    np.testing.assert_almost_equal(F3E, F3e)

    return 1


def test_Form_No8_edge_forms():
    if rAnk == mAster_rank:
        print("EEE [test_Form_No8_edge_forms] ...... ", flush=True)

    return 1


def test_Form_No9_node_forms():
    if rAnk == mAster_rank:
        print("NNN [test_Form_No9_node_forms] ...... ", flush=True)

    return 1




def test_Form_No10_standard_form_dofs():
    if rAnk == mAster_rank:
        print("DOF [test_Form_No10_standard_form_dofs] ...... ", flush=True)


    if rAnk == mAster_rank:
        load = random.randint(100,499)
    else:
        load = None
    load = cOmm.bcast(load, root=mAster_rank)
    FC = random_3D_FormCaller_of_total_load_around(load)

    f0 = FC('0-f', is_hybrid=False)
    dofs = f0.dofs
    for i in dofs:
        assert i in dofs, f"must be the case"
        D = dofs[i]
        local_positions = D.positions
        for E_I in local_positions:
            E, I = E_I
            assert f0.numbering.gathering[E][I] == i, f"must be the case."
    return 1


def test_Form_No11_reconstruction_matrices():
    """"""
    if rAnk == mAster_rank:
        load = random.randint(100,500)
        IH = [True, False][random.randint(0,1)]
        print(f"~~~ [test_Form_No11_reconstruction_matrices] @ load= {load}... ", flush=True)
    else:
        load = None
        IH = None
    load, IH = cOmm.bcast([load, IH], root=mAster_rank)
    FC = random_3D_FormCaller_of_total_load_around(load)
    a = random.random()
    b = random.random()
    c = random.random()
    def P(t, x, y, z): return - 1.76*np.pi * np.sin(2.1*np.pi*x) * np.cos(1.11*np.pi*y) * np.sin(a * 2.12*np.pi*z) + t/2
    def Q(t, x, y, z): return np.pi * np.cos(b * 3.52*np.pi*x) * np.sin(2*np.pi*y) * np.cos(1.31*np.pi*z)+ t
    def R(t, x, y, z): return 3.15 * np.pi * np.cos(4.52*np.pi*x) * np.sin(1.5*np.pi*y) + t + np.cos(c * 2.11*np.pi*z)

    scalar = FC('scalar', P)
    vector = FC('vector', (P, Q, R))
    xi = np.linspace(-1,1,7)
    et = np.linspace(-1,1,6)
    sg = np.linspace(-1,1,5)
    #------------- 1-form ------------------------------
    f = FC('1-f', is_hybrid=IH)
    f.TW.func.do.set_func_body_as(vector)
    f.TW.current_time = 1
    f.TW.do.push_all_to_instant()
    f.discretize()

    R = f.reconstruct(xi, et, sg, ravel=True)[1]
    MR = f.do.make_reconstruction_matrix_on_grid(xi, et, sg)
    for i in MR:
        U, V, W = MR[i] @ f.cochain.local[i]
        u, v, w = R[i]
        np.testing.assert_almost_equal(np.max(np.abs(U- u)), 0)
        np.testing.assert_almost_equal(np.max(np.abs(v- V)), 0)
        np.testing.assert_almost_equal(np.max(np.abs(w- W)), 0)

    #------------- 2-form ------------------------------
    f = FC('2-f', is_hybrid=IH)
    f.TW.func.do.set_func_body_as(vector)
    f.TW.current_time = 2
    f.TW.do.push_all_to_instant()
    f.discretize()

    R = f.reconstruct(xi, et, sg, ravel=True)[1]
    MR = f.do.make_reconstruction_matrix_on_grid(xi, et, sg)
    for i in MR:
        U, V, W = MR[i] @ f.cochain.local[i]
        u, v, w = R[i]
        np.testing.assert_almost_equal(np.max(np.abs(U- u)), 0)
        np.testing.assert_almost_equal(np.max(np.abs(v- V)), 0)
        np.testing.assert_almost_equal(np.max(np.abs(w- W)), 0)

    #------------ 0-form ----------------------------------------------------
    f = FC('0-f', is_hybrid=IH)
    f.TW.func.do.set_func_body_as(scalar)
    f.TW.current_time = 3
    f.TW.do.push_all_to_instant()
    f.discretize()
    R = f.reconstruct(xi, et, sg, ravel=True)[1]
    MR = f.do.make_reconstruction_matrix_on_grid(xi, et, sg)
    for i in MR:
        r = MR[i] @ f.cochain.local[i]
        np.testing.assert_almost_equal(np.max(np.abs(r - R[i][0])), 0)

    #------------ 3-form -------------------------------------------------------
    f = FC('3-f', is_hybrid=IH)
    f.TW.func.do.set_func_body_as(scalar)
    f.TW.current_time = 4
    f.TW.do.push_all_to_instant()
    f.discretize()
    R = f.reconstruct(xi, et, sg, ravel=True)[1]
    MR = f.do.make_reconstruction_matrix_on_grid(xi, et, sg)
    for i in MR:
        r = MR[i] @ f.cochain.local[i]
        np.testing.assert_almost_equal(np.max(np.abs(r - R[i][0])), 0)

    return 1


if __name__ == '__main__':
    # mpiexec -n 4 python _3dCSCG\tests\unittests\standard_forms.py
    test_Form_NOx_cross_product_3()