# -*- coding: utf-8 -*-
"""
For standard forms only.

"""
import sys

if './' not in sys.path: sys.path.append('./')

import random

from root.config.main import *

from scipy.sparse import linalg as spspalinalg

from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller, ExactSolutionSelector
from components.exceptions import ThreeDimensionalTransfiniteInterpolationError
from tests.objects.CSCG._3d.randObj.form_caller import random_mesh_and_space_of_total_load_around
from tests.objects.CSCG._3d.randObj.form_caller import random_FormCaller_of_total_load_around
from tests.objects.CSCG._3d.randObj.field import random_vector

from tools.miLinearAlgebra.linearSystem.main import LinearSystem
from tools.elementwiseCache.dataStructures.operators.bmat.main import bmat
from tools.elementwiseCache.dataStructures.operators.concatenate.main import concatenate

def test_Form_NO1_discretization_and_reconstruction():
    """"""
    if RANK == MASTER_RANK:
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
        f0 = FC('0-f', hybrid=False)
        f1 = FC('1-f', hybrid=False)
        f2 = FC('2-f', hybrid=False)
        f3 = FC('3-f', hybrid=False)
        f0.CF = scalar
        f0.CF.current_time = 0
        f0.discretize()
        assert f0.error.L() < 0.00003
        f1.CF = vector
        f1.CF.current_time = 0
        f1.discretize()
        assert f1.error.L() < 0.0004
        f2.CF = vector
        f2.CF.current_time = 0
        f2.discretize()
        assert f2.error.L() < 0.0005
        f3.CF = scalar
        f3.CF.current_time = 0
        f3.discretize()
        assert f3.error.L() < 0.0004
    except ThreeDimensionalTransfiniteInterpolationError:
        if RANK == MASTER_RANK:
            print("    skip LDC mesh ...... ", flush=True)
        pass

    mesh = MeshGenerator('crazy_periodic', c=0.1)([3, 4, 5], EDM='debug')
    space = SpaceInvoker('polynomials')([('Lobatto', 5), ('Lobatto', 4), ('Lobatto', 3)])
    FC = FormCaller(mesh, space)

    scalar = FC('scalar', p)
    vector = FC('vector', (u,v,w))

    f0 = FC('0-f', hybrid=False)
    f1 = FC('1-f', hybrid=False)
    f2 = FC('2-f', hybrid=False)
    f3 = FC('3-f', hybrid=False)

    f0.CF = scalar
    f0.CF.current_time = 0
    f0.discretize()
    assert f0.error.L() < 0.0004
    f1.CF = vector
    f1.CF.current_time = 0
    f1.discretize()
    assert f1.error.L() < 0.004
    f2.CF = vector
    f2.CF.current_time = 0
    f2.discretize()
    assert f2.error.L() < 0.0041
    f3.CF = scalar
    f3.CF.current_time = 0
    f3.discretize()
    assert f3.error.L() < 0.003

    f0.cochain.globe = f0.cochain.globe
    f1.cochain.globe = f1.cochain.globe
    f2.cochain.globe = f2.cochain.globe
    f3.cochain.globe = f3.cochain.globe
    f0.CF = scalar
    f0.CF.current_time = 0
    f0.discretize()
    assert f0.error.L() < 0.0004
    f1.CF = vector
    f1.CF.current_time = 0
    f1.discretize()
    assert f1.error.L() < 0.004
    f2.CF = vector
    f2.CF.current_time = 0
    f2.discretize()
    assert f2.error.L() < 0.0041
    f3.CF = scalar
    f3.CF.current_time = 0
    f3.discretize()
    assert f3.error.L() < 0.003

    f0 = FC('0-f', hybrid=True)
    f1 = FC('1-f', hybrid=True)
    f2 = FC('2-f', hybrid=True)
    f3 = FC('3-f', hybrid=True)

    f0.CF = scalar
    f0.CF.current_time = 0
    f0.discretize()
    assert f0.error.L() < 0.0004
    f1.CF = vector
    f1.CF.current_time = 0
    f1.discretize()
    assert f1.error.L() < 0.004
    f2.CF = vector
    f2.CF.current_time = 0
    f2.discretize()
    assert f2.error.L() < 0.0041
    f3.CF = scalar
    f3.CF.current_time = 0
    f3.discretize()
    assert f3.error.L() < 0.003

    f0.cochain.globe = f0.cochain.globe
    f1.cochain.globe = f1.cochain.globe
    f2.cochain.globe = f2.cochain.globe
    f3.cochain.globe = f3.cochain.globe
    f0.CF = scalar
    f0.CF.current_time = 0
    f0.discretize()
    assert f0.error.L() < 0.0004
    f1.CF = vector
    f1.CF.current_time = 0
    f1.discretize()
    assert f1.error.L() < 0.004
    f2.CF = vector
    f2.CF.current_time = 0
    f2.discretize()
    assert f2.error.L() < 0.0041
    f3.CF = scalar
    f3.CF.current_time = 0
    f3.discretize()
    assert f3.error.L() < 0.003


    mesh = MeshGenerator('bridge_arch_cracked',)([3,2,4], EDM='debug')
    space = SpaceInvoker('polynomials')([('Lobatto', 3), ('Lobatto', 4), ('Lobatto', 2)])
    FC = FormCaller(mesh, space)

    scalar = FC('scalar', p)
    vector = FC('vector', (u,v,w))
    f0 = FC('0-f', hybrid=False)
    f1 = FC('1-f', hybrid=False)
    f2 = FC('2-f', hybrid=False)
    f3 = FC('3-f', hybrid=False)

    f0.CF = scalar
    f0.CF.current_time = 0
    f0.discretize()
    assert f0.error.L() < 0.01
    f1.CF = vector
    f1.CF.current_time = 0
    f1.discretize()
    assert f1.error.L() < 0.05
    f2.CF = vector
    f2.CF.current_time = 0
    f2.discretize()
    assert f2.error.L() < 0.07
    f3.CF = scalar
    f3.CF.current_time = 0
    f3.discretize()
    assert f3.error.L() < 0.05


    mesh = MeshGenerator('crazy', c=0.1)([5,5,5], EDM=None)
    space = SpaceInvoker('polynomials')([('Lobatto', 4), ('Lobatto', 4), ('Lobatto', 4)])
    FC = FormCaller(mesh, space)

    scalar = FC('scalar', p)
    vector = FC('vector', (u,v,w))
    f0 = FC('0-f', hybrid=False)
    f1 = FC('1-f', hybrid=False)
    f2 = FC('2-f', hybrid=False)
    f3 = FC('3-f', hybrid=False)

    f0.CF = scalar
    f0.CF.current_time = 0
    f0.discretize()
    assert f0.error.L() < 0.00005
    f1.CF = vector
    f1.CF.current_time = 0
    f1.discretize()
    assert f1.error.L() < 0.0008
    f2.CF = vector
    f2.CF.current_time = 0
    f2.discretize()
    assert f2.error.L() < 0.0008
    f3.CF = scalar
    f3.CF.current_time = 0
    f3.discretize()
    assert f3.error.L() < 0.0006

    return 1


def test_Form_NO1a_discretization_and_reconstruction():
    """"""
    if RANK == MASTER_RANK:
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
        f0 = FC('0-f', hybrid=False)
        f1 = FC('1-f', hybrid=False)
        f2 = FC('2-f', hybrid=False)
        f3 = FC('3-f', hybrid=False)
        f0.CF = scalar
        f0.CF.current_time = 0
        f0.discretize()
        assert f0.error.L() < 0.00003
        f1.CF = vector
        f1.CF.current_time = 0
        f1.discretize()
        assert f1.error.L() < 0.0004
        f2.CF = vector
        f2.CF.current_time = 0
        f2.discretize()
        assert f2.error.L() < 0.0005
        f3.CF = scalar
        f3.CF.current_time = 0
        f3.discretize()
        assert f3.error.L() < 0.0004
    except ThreeDimensionalTransfiniteInterpolationError:
        if RANK == MASTER_RANK:
            print("    skip LDC mesh ...... ", flush=True)
        pass

    mesh = MeshGenerator('crazy_periodic', c=0.1)([3, 4, 5], EDM='chaotic')
    space = SpaceInvoker('polynomials')([('Lobatto', 5), ('Lobatto', 4), ('Lobatto', 3)])
    FC = FormCaller(mesh, space)

    scalar = FC('scalar', p)
    vector = FC('vector', (u,v,w))

    f0 = FC('0-f', hybrid=False)
    f1 = FC('1-f', hybrid=False)
    f2 = FC('2-f', hybrid=False)
    f3 = FC('3-f', hybrid=False)

    f0.CF = scalar
    f0.CF.current_time = 0
    f0.discretize()
    assert f0.error.L() < 0.0004
    f1.CF = vector
    f1.CF.current_time = 0
    f1.discretize()
    assert f1.error.L() < 0.004
    f2.CF = vector
    f2.CF.current_time = 0
    f2.discretize()
    assert f2.error.L() < 0.0041
    f3.CF = scalar
    f3.CF.current_time = 0
    f3.discretize()
    assert f3.error.L() < 0.003

    f0.cochain.globe = f0.cochain.globe
    f1.cochain.globe = f1.cochain.globe
    f2.cochain.globe = f2.cochain.globe
    f3.cochain.globe = f3.cochain.globe
    f0.CF = scalar
    f0.CF.current_time = 0
    f0.discretize()
    assert f0.error.L() < 0.0004
    f1.CF= vector
    f1.CF.current_time = 0
    f1.discretize()
    assert f1.error.L() < 0.004
    f2.CF = vector
    f2.CF.current_time = 0
    f2.discretize()
    assert f2.error.L() < 0.0041
    f3.CF = scalar
    f3.CF.current_time = 0
    f3.discretize()
    assert f3.error.L() < 0.003

    f0 = FC('0-f', hybrid=True)
    f1 = FC('1-f', hybrid=True)
    f2 = FC('2-f', hybrid=True)
    f3 = FC('3-f', hybrid=True)

    f0.CF = scalar
    f0.CF.current_time = 0
    f0.discretize()
    assert f0.error.L() < 0.0004
    f1.CF = vector
    f1.CF.current_time = 0
    f1.discretize()
    assert f1.error.L() < 0.004
    f2.CF = vector
    f2.CF.current_time = 0
    f2.discretize()
    assert f2.error.L() < 0.0041
    f3.CF = scalar
    f3.CF.current_time = 0
    f3.discretize()
    assert f3.error.L() < 0.003

    f0.cochain.globe = f0.cochain.globe
    f1.cochain.globe = f1.cochain.globe
    f2.cochain.globe = f2.cochain.globe
    f3.cochain.globe = f3.cochain.globe
    f0.CF = scalar
    f0.CF.current_time = 0
    f0.discretize()
    assert f0.error.L() < 0.0004
    f1.CF = vector
    f1.CF.current_time = 0
    f1.discretize()
    assert f1.error.L() < 0.004
    f2.CF = vector
    f2.CF.current_time = 0
    f2.discretize()
    assert f2.error.L() < 0.0041
    f3.CF = scalar
    f3.CF.current_time = 0
    f3.discretize()
    assert f3.error.L() < 0.003


    mesh = MeshGenerator('bridge_arch_cracked',)([3,2,4], EDM='chaotic')
    space = SpaceInvoker('polynomials')([('Lobatto', 3), ('Lobatto', 4), ('Lobatto', 2)])
    FC = FormCaller(mesh, space)

    scalar = FC('scalar', p)
    vector = FC('vector', (u,v,w))
    f0 = FC('0-f', hybrid=False)
    f1 = FC('1-f', hybrid=False)
    f2 = FC('2-f', hybrid=False)
    f3 = FC('3-f', hybrid=False)

    f0.CF = scalar
    f0.CF.current_time = 0
    f0.discretize()
    assert f0.error.L() < 0.01
    f1.CF = vector
    f1.CF.current_time = 0
    f1.discretize()
    assert f1.error.L() < 0.05
    f2.CF = vector
    f2.CF.current_time = 0
    f2.discretize()
    assert f2.error.L() < 0.07
    f3.CF = scalar
    f3.CF.current_time = 0
    f3.discretize()
    assert f3.error.L() < 0.05

    return 1


def test_Form_NO1b_trace_form_Rd_and_Rc():
    """"""
    if RANK == MASTER_RANK:
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
    t0.CF = flux
    t0.CF.current_time = t
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
    t2.CF = flux
    t2.CF.current_time = t
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

    if RANK == MASTER_RANK:
        load = random.randint(500,1000)
        t = random.random()
    else:
        load, t = None, None
    load, t = COMM.bcast([load, t], root=MASTER_RANK)
    mesh, space = random_mesh_and_space_of_total_load_around(load, exclude_periodic = True)

    FC = FormCaller(mesh, space)
    flux = FC('scalar', p)
    velo = FC('vector', (uuu, vvv, www))

    t0 = FC('0-t')
    t1 = FC('1-t')
    t2 = FC('2-t')

    f0 = FC('0-f', hybrid=True)
    f1 = FC('1-f', hybrid=True)
    f2 = FC('2-f', hybrid=True)

    S0 = t0.matrices.selective
    S1 = t1.matrices.selective
    S2 = t2.matrices.selective

    # t0 & f0, discretization and selective matrix
    t0.CF = flux
    t0.CF.current_time = t
    t0.discretize()
    f0.CF = flux
    f0.CF.current_time = t
    f0.discretize()
    t0_local = t0.cochain.local
    f0_local = f0.cochain.local
    for i in t0_local:
        A0 = t0_local[i] -  S0[i] @ f0_local[i]
        np.testing.assert_array_almost_equal(A0, 0)

    # t1 & f1, discretization and selective matrix
    t1.CF = velo
    t1.CF.current_time = t
    t1.discretize()
    f1.CF = velo
    f1.CF.current_time = t
    f1.discretize()
    t1_local = t1.cochain.local
    f1_local = f1.cochain.local
    for i in t1_local:
        A1 = t1_local[i] -  S1[i] @ f1_local[i]
        np.testing.assert_array_almost_equal(A1, 0)

    # t2 & f2, discretization and selective matrix
    t2.CF = velo
    t2.CF.current_time = t
    t2.discretize()
    f2.CF = velo
    f2.CF.current_time = t
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
    if RANK == MASTER_RANK:
        print("*** [test_Form_NO2_mass_matrix] ...... ", flush=True)

    mesh = MeshGenerator('crazy', c=0.24)([2, 2, 2], EDM='debug')
    space = SpaceInvoker('polynomials')([('Lobatto', 2), ('Lobatto', 2), ('Lobatto', 2)])
    FC = FormCaller(mesh, space)

    f0 = FC('0-f', hybrid=False)
    f1 = FC('1-f', hybrid=False)
    f2 = FC('2-f', hybrid=True)
    f3 = FC('3-f', hybrid=True)

    benchmark0 = np.array([
        3.54729207e-04,  1.95995397e-04, -9.05682361e-05,  1.95995397e-04,
        1.13207590e-04, -5.03123358e-05, -9.05682361e-05, -5.03123358e-05,
        2.20104383e-05,  1.95995397e-04,  1.13207590e-04, -5.03123358e-05,
        1.13207590e-04,  6.90036501e-05, -2.91740562e-05, -5.03123358e-05,
       -2.91740562e-05,  1.19231437e-05, -9.05682361e-05, -5.03123358e-05,
        2.20104383e-05, -5.03123358e-05, -2.91740562e-05,  1.19231437e-05,
        2.20104383e-05,  1.19231437e-05, -4.62962963e-06])
    benchmark1 = np.array([
        2.45827370e-02, -2.27848002e-03, 1.33590655e-02, -9.95744998e-04,
        -5.99257863e-03, 5.78395986e-04, 1.33590655e-02, -9.95744998e-04,
        7.47441051e-03, -3.20884734e-04, -3.09020211e-03, 3.18213891e-04,
        -5.99257863e-03, 5.78395986e-04, -3.09020211e-03, 3.18213891e-04,
        1.29851550e-03, -1.69522326e-04, -1.96278457e-03, -2.93898015e-03,
        1.14722269e-04, 2.85675816e-04, -1.37729255e-04, 9.13028259e-05,
        -1.70509316e-03, -2.83076200e-03, 1.11087334e-04, 2.55836077e-04,
        1.88421423e-04, 8.84774860e-05, 8.00331409e-04, 1.21016986e-03,
        -2.13948002e-05, -1.00747199e-04, 1.54140389e-04, -5.14863383e-05,
        -1.96278457e-03, -2.93898015e-03, 1.14722269e-04, -1.70509316e-03,
        -2.83076200e-03, 1.11087334e-04, 8.00331409e-04, 1.21016986e-03,
        -2.13948002e-05, 2.85675816e-04, -1.37729255e-04, 9.13028259e-05,
        2.55836077e-04, 1.88421423e-04, 8.84774860e-05, -1.00747199e-04,
        1.54140389e-04, -5.14863383e-05])
    benchmark2 = np.array([
        1.58529903e+00, 7.75026284e-01, -3.26457821e-01, -1.70723220e-01,
        -9.37333827e-02, 5.23050991e-02, -1.70723220e-01, -9.37333827e-02,
        5.23050991e-02, 3.40535298e-02, 8.09800992e-03, -5.94674600e-03,
        2.16808181e-01, -3.01069942e-02, 3.57773956e-01, -9.06827115e-03,
        -1.14611365e-02, -9.53669458e-03, 2.54466757e-02, -8.20605629e-04,
        3.73213896e-02, 1.08921597e-02, 2.90710672e-03, -3.26394743e-03,
        2.16808181e-01, -3.01069942e-02, 2.54466757e-02, -8.20605629e-04,
        3.57773956e-01, -9.06827115e-03, 3.73213896e-02, 1.08921597e-02,
        -1.14611365e-02, -9.53669458e-03, 2.90710672e-03, -3.26394743e-03])
    benchmark3 = np.array([
        78.35846174, -17.15884714, -17.15884714, 1.02217158,
         -17.15884714, 1.02217158, 1.02217158, -0.86973663])

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
    if RANK == MASTER_RANK:
        print("*** [test_Form_NO3_incidence_matrices] ...... ", flush=True)

    mesh = MeshGenerator('crazy', c=0.0)([2, 3, 4], EDM='debug')
    space = SpaceInvoker('polynomials')([('Lobatto', 8), ('Lobatto', 7), ('Lobatto', 6)])
    FC = FormCaller(mesh, space)
    es = ExactSolutionSelector(mesh)('icpsNS:sincosRD')
    np.testing.assert_almost_equal(es.helicity(0), -12.566370614359162)
    np.testing.assert_almost_equal(es.kinetic_energy(0), 1.4999999999999998)
    np.testing.assert_almost_equal(es.enstrophy(0), 59.217626406536155)
    # so kinetic_energy = 0.5 * (L^2 norm of velocity)**2
    np.testing.assert_almost_equal(
        es.___Pr_compute_Ln_norm_of___('velocity', time=0, n=2), np.sqrt(3))

    f0 = FC('0-f', hybrid=True)

    f0.CF = es.pressure
    f0.CF.current_time = 0
    f0.discretize()
    assert f0.error.L() < 5.8e-7

    f1 = f0.coboundary()
    f1.CF = es.gradient_of_pressure
    f1.CF.current_time = 0
    assert f1.error.L() < 3.6e-5

    f1.CF = es.velocity
    f1.CF.current_time = 0
    f1.discretize()
    f2 = f1.coboundary()
    f2.CF = es.vorticity
    f2.CF.current_time = 0
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

    f3.CF = es.divergence_of_velocity
    f3.CF.current_time = 0
    L_inf = f3.error.L(n='infinity')
    assert L_inf < 1e-8

    # we test curl of vorticity for incompressible NS ----------- BELOW ------------------------------------------------
    f1 = FC('1-f', hybrid=False)
    f1.CF = es.vorticity
    f1.CF.current_time = 0
    f1.discretize()
    f2 = f1.coboundary()
    f2.CF = es.curl_of_vorticity
    f2.CF.current_time = 0
    assert f2.error.L() < 0.00045

    mesh = MeshGenerator('crazy', c=0.0)([6, 6, 6], EDM='debug')
    space = SpaceInvoker('polynomials')([('Lobatto', 8), ('Lobatto', 8), ('Lobatto', 8)])
    FC = FormCaller(mesh, space)
    es = ExactSolutionSelector(mesh)('icpsNS:sincosRD')

    f1 = FC('1-f', hybrid=True)
    f1.CF = es.vorticity
    f1.CF.current_time = 0
    f1.discretize()
    f2 = f1.coboundary()
    f2.CF = es.curl_of_vorticity
    f2.CF.current_time = 0
    assert f2.error.L() < 1e-7

    # we test numerical gradient, curl and div with coboundary of standard forms------------------------------------
    if RANK == MASTER_RANK:
        c = random.random() / 10
        if c < 0.05:
            c = 0
        else:
            pass
    else:
        c = None
    c = COMM.bcast(c, root=MASTER_RANK)
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

    f0 = FC('0-f', hybrid=False)
    f0.CF= S
    f0.CF.current_time = t
    f0.discretize()
    f1 = f0.coboundary()
    f1.CF = gS
    f1.CF.current_time = t
    assert f1.error.L() < 0.0065

    f1 = FC('1-f', hybrid=False)
    f1.CF = V
    f1.CF.current_time = t
    f1.discretize()
    f2 = f1.coboundary()
    f2.CF = cV
    f2.CF.current_time = t
    assert f2.error.L() < 0.002

    f2 = FC('2-f', hybrid=False)
    f2.CF = V
    f2.CF.current_time = t
    f2.discretize()
    f3 = f2.coboundary()
    f3.CF = dV
    f3.CF.current_time = t
    assert f3.error.L() < 0.0025

    # +++++++++++++++++++++++++++++++++++++++++ ABOVE ++++++++++++++++++++++++++++++++++++++++++++++

    return 1


def test_Form_NO4_cross_product_1():
    """
    Unittests for the special method ``cross_product`` of the standard 1-form.
    """
    if RANK == MASTER_RANK:
        print("*** [test_Form_NO4_cross_product_1] ...... ", flush=True)

    mesh = MeshGenerator('crazy', c=0.0)([3, 3, 4], EDM='debug')
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

    w.CF = vw
    w.CF.current_time = 0
    w.discretize()

    # represent (w times u, e) or (w^1 wedge u^1, e^1)
    MW = w.special.cross_product_1f__ip_1f(u, x)
    M = x.matrices.mass

    u.CF = vu
    u.CF.current_time = 0
    u.discretize()
    _xcl = {}
    for i in mesh.elements:
        P = spspalinalg.inv(M[i].tocsc()) @ MW[i]
        _xcl[i] = P @ u.cochain.local[i]
    x.cochain.local = _xcl

    x.CF = vx
    x.CF.current_time = 0
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

    w.CF = vw
    w.CF.current_time = 0
    w.discretize()

    # represent (w times u, e) or (w^1 wedge u^1, e^1)
    MW = w.special.cross_product_1f__ip_1f(u, x)
    M = x.matrices.mass

    u.CF = vu
    u.CF.current_time = 0
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

    x.CF = vx
    x.CF.current_time = 0
    assert x.error.L() < 0.025, f"{x.error.L()}"

    return 1


def test_Form_NO5_cross_product_2():
    """
    Unittests for the special method ``cross_product`` of the standard 2-form.
    """
    if RANK == MASTER_RANK:
        print("*** [test_Form_NO5_cross_product_2] ...... ", flush=True)

    mesh = MeshGenerator('crazy', c=0.0, bounds=([0,1],[0,1],[0,1]))([4,3,5], EDM='debug')
    space = SpaceInvoker('polynomials')([('Lobatto', 3), ('Lobatto', 4), ('Lobatto', 2)])
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

    w.CF = vw
    w.CF.current_time = 0
    w.discretize()

    # represent (w times u, e) or (star w^2 wedge star u^2, e^2)
    MW = w.special.cross_product_2f__ip_2f(u, x)
    M = x.matrices.mass

    u.CF = vu
    u.CF.current_time = 0
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

    x.CF = vx
    x.CF.current_time = 0
    assert x.error.L() < 0.015

    return 1


def test_Form_NOx_cross_product_3():
    """
    Unittests for the special method ``cross_product`` of the standard 1-form.
    """
    if RANK == MASTER_RANK:
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

    w.CF = vw
    w.CF.current_time = 0
    w.discretize()

    # represent (w1 times u2, e2)
    MW = w.special.cross_product_2f__ip_2f(u, x)
    M = x.matrices.mass

    u.CF = vu
    u.CF.current_time = 0
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

    x.CF = vx
    x.CF.current_time = 0
    assert x.error.L() < 0.0055

    return 1


def test_Form_NOx1_cross_product_4():
    """Unittests for the special method ``curl_self_cross_product_self__ip_2f`` of the standard 1-form.
    """
    if RANK == MASTER_RANK:
        print("*** [test_Form_NOx1_cross_product_4] ...... ", flush=True)

    mesh = MeshGenerator('cuboid', region_layout=[2,2,2])([5,5,5])
    space = SpaceInvoker('polynomials')([2,2,2])
    FC = FormCaller(mesh, space)

    def u0(t, x, y, z): return np.sin(0.87*np.pi*x) * np.sin(2*y) * np.cos(0.55*np.pi*z) + t
    def u1(t, x, y, z): return np.cos(2.68*x) * np.sin(0.75*np.pi*y) * np.sin(0.82*np.pi*z) + t
    def u2(t, x, y, z): return np.sin(0.854*np.pi*x) * np.cos(0.56*np.pi*y) * np.sin(1.987*z) + t

    u = FC('vector', func=(u0, u1, u2))
    w = u.numerical.curl
    wXu = w.do.cross_product(u)

    U = FC('1-f')
    E = FC('2-f')
    U.CF = u
    U.CF.current_time = 0
    U.discretize()

    CPV = U.special.curl_self_cross_product_self__ip_2f(E)
    M = E.matrices.mass

    if mesh.elements.num > 0:
        i = mesh.elements.indices[0]
        invM = spspalinalg.inv(M[i].tocsc())
    else:
        invM = None

    Ecl = dict()
    E21 = U.matrices.incidence

    for i in mesh.elements:
        v = CPV[i].toarray()[:,0]
        Ecl[i] = invM @ v

        result = np.sum((E21[i] @ U.cochain.local[i]) * v)
        np.testing.assert_almost_equal(result, 0)

    E.cochain.local = Ecl
    E.CF = wXu
    E.CF.current_time = 0
    assert E.error.L() < 0.01

    return 1


def test_Form_NOx2_F_dot_G_times_H():
    """"""
    if RANK == MASTER_RANK:
        print("*** [test_Form_NOx2_F_dot_G_times_H] ...... ", flush=True)

    mesh = MeshGenerator('cuboid', region_layout=[2,2,2])([5,5,5])
    space = SpaceInvoker('polynomials')([2,2,2])
    FC = FormCaller(mesh, space)
    uV = random_vector(mesh)
    BV = random_vector(mesh)
    jV = random_vector(mesh)
    u = FC('2-f')
    j = FC('1-f')
    B = FC('2-f')

    uXB = uV.do.cross_product(BV)
    uXB_dot_j = uXB.do.inner_product(jV)

    jXB = jV.do.cross_product(BV)
    jXB_dot_u = jXB.do.inner_product(uV)

    uXB_dot_j.current_time = 0
    jXB_dot_u.current_time = 0

    xi = np.random.rand(3,3,3)
    et = np.random.rand(3,3,3)
    sg = np.random.rand(3,3,3)

    N0 = uXB_dot_j.do.compute_Ln_norm(n=1)
    N1 = jXB_dot_u.do.compute_Ln_norm(n=1)
    np.testing.assert_almost_equal(N0 + N1, 0)

    _, Val0 = uXB_dot_j.do.reconstruct(xi, et, sg)
    _, Val1 = jXB_dot_u.do.reconstruct(xi, et, sg)
    for i in Val0:
        np.testing.assert_array_almost_equal(Val0[i][0] + Val1[i][0], 0)

    u.CF = uV
    u.CF.current_time = 0
    u.discretize()
    B.CF = BV
    B.CF.current_time = 0
    B.discretize()
    j.CF = jV
    j.CF.current_time = 0
    j.discretize()

    CM0 = u.special.cross_product_2f__ip_1f(B, j)
    CM1 = j.special.cross_product_2f__ip_2f(B, u)

    for i in CM0:
        val0 = np.einsum('ij, i, j ->', CM0[i].toarray(), j.cochain.local[i], B.cochain.local[i], optimize='optimal')
        val1 = np.einsum('ij, i, j ->', CM1[i].toarray(), u.cochain.local[i], B.cochain.local[i], optimize='optimal')
        np.testing.assert_almost_equal(val0 + val1, 0)

    return 1



def test_Form_No7_with_other_element_numbering_AUTO():
    if RANK == MASTER_RANK:
        print("*** [test_Form_No7_with_other_element_numbering---AUTO] ...... ", flush=True)

    def u(t, x, y, z): return np.cos(np.pi*x) + np.sin(np.pi*y) * np.sin(np.pi*z-0.125)**2 + t/2
    def v(t, x, y, z): return np.sin(np.pi*x) + np.sin(np.pi*y) * np.sin(np.pi*z-0.125)**2 + t/2
    def w(t, x, y, z): return np.sin(np.pi*x) + np.cos(np.pi*y) * np.cos(np.pi*z-0.125)**2 + t
    def p(t, x, y, z): return np.cos(np.pi*x) + np.sin(np.pi*y) * np.sin(np.pi*z-0.125)**2 + t/2

    if RANK == MASTER_RANK:
        t = random.uniform(-100,100)
    else:
        t = None
    t = COMM.bcast(t, root=MASTER_RANK)

    if RANK == MASTER_RANK:
        i = random.randint(2,4)
        j = random.randint(2,4)
        k = random.randint(2,4)
    else:
        i, j, k = None, None, None
    i, j, k = COMM.bcast([i, j, k], root=MASTER_RANK)
    space = SpaceInvoker('polynomials')([('Lobatto', i), ('Lobatto', j), ('Lobatto', k)])

    if RANK == MASTER_RANK:
        i = random.randint(5,8)
        j = random.randint(5,9)
        k = random.randint(6,7)
    else:
        i, j, k = None, None, None
    i, j, k = COMM.bcast([i, j, k], root=MASTER_RANK)
    MESH = MeshGenerator('bridge_arch_cracked')([i, j, k], EDM='debug', show_info=False)
    mesh = MeshGenerator('bridge_arch_cracked')([i, j, k], EDM=None, show_info=False)
    MQ = MESH.___PRIVATE_element_division_and_numbering_quality___()[0]
    mQ = mesh.___PRIVATE_element_division_and_numbering_quality___()[0]
    assert MQ <= mQ, "We expect this!"

    FC = FormCaller(mesh, space)
    scalar = FC('scalar', p)
    vector = FC('vector', (u,v,w))
    f0 = FC('0-f', hybrid=False)
    f1 = FC('1-f', hybrid=False)
    f2 = FC('2-f', hybrid=False)
    f3 = FC('3-f', hybrid=False)
    f0.CF = scalar
    f0.CF.current_time = t
    f0.discretize()
    F0E = f0.error.L()
    f1.CF = vector
    f1.CF.current_time = t
    f1.discretize()
    F1E = f1.error.L()
    f2.CF = vector
    f2.CF.current_time = t
    f2.discretize()
    F2E = f2.error.L()
    f3.CF = scalar
    f3.CF.current_time = t
    f3.discretize()
    F3E = f3.error.L()


    FC = FormCaller(MESH, space)
    scalar = FC('scalar', p)
    vector = FC('vector', (u,v,w))
    f0 = FC('0-f', hybrid=False)
    f1 = FC('1-f', hybrid=False)
    f2 = FC('2-f', hybrid=False)
    f3 = FC('3-f', hybrid=False)
    f0.CF = scalar
    f0.CF.current_time = t
    f0.discretize()
    F0e = f0.error.L()
    f1.CF = vector
    f1.CF.current_time = t
    f1.discretize()
    F1e = f1.error.L()
    f2.CF = vector
    f2.CF.current_time = t
    f2.discretize()
    F2e = f2.error.L()
    f3.CF = scalar
    f3.CF.current_time = t
    f3.discretize()
    F3e = f3.error.L()

    np.testing.assert_almost_equal(F0E, F0e)
    np.testing.assert_almost_equal(F1E, F1e)
    np.testing.assert_almost_equal(F2E, F2e)
    np.testing.assert_almost_equal(F3E, F3e)

    return 1


def test_Form_No10_standard_form_dofs():
    if RANK == MASTER_RANK:
        print("DOF [test_Form_No10_standard_form_dofs] ...... ", flush=True)


    if RANK == MASTER_RANK:
        load = random.randint(100,499)
    else:
        load = None
    load = COMM.bcast(load, root=MASTER_RANK)
    FC = random_FormCaller_of_total_load_around(load)

    f0 = FC('0-f', hybrid=False)
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
    if RANK == MASTER_RANK:
        load = random.randint(100,200)
        IH = [True, False][random.randint(0,1)]
        print(f"~~~ [test_Form_No11_reconstruction_matrices] @ load= {load}... ", flush=True)
    else:
        load = None
        IH = None
    load, IH = COMM.bcast([load, IH], root=MASTER_RANK)
    FC = random_FormCaller_of_total_load_around(load, mesh_pool=('crazy',))
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
    f = FC('1-f', hybrid=IH)
    f.CF = vector
    f.CF.current_time = 1
    f.discretize()

    R = f.reconstruct(xi, et, sg, ravel=True)[1]
    MR = f.do.make_reconstruction_matrix_on_grid(xi, et, sg)
    for i in MR:
        U, V, W = MR[i] @ f.cochain.local[i]
        u, v, w = R[i]
        np.testing.assert_almost_equal(np.max(np.abs(U - u)), 0)
        np.testing.assert_almost_equal(np.max(np.abs(v - V)), 0)
        np.testing.assert_almost_equal(np.max(np.abs(w - W)), 0)

    #------------- 2-form ------------------------------
    f = FC('2-f', hybrid=IH)
    f.CF = vector
    f.CF.current_time = 2
    f.discretize()

    R = f.reconstruct(xi, et, sg, ravel=True)[1]
    MR = f.do.make_reconstruction_matrix_on_grid(xi, et, sg)
    for i in MR:
        U, V, W = MR[i] @ f.cochain.local[i]
        u, v, w = R[i]
        np.testing.assert_almost_equal(np.max(np.abs(U - u)), 0)
        np.testing.assert_almost_equal(np.max(np.abs(v - V)), 0)
        np.testing.assert_almost_equal(np.max(np.abs(w - W)), 0)

    #------------ 0-form ----------------------------------------------------
    f = FC('0-f', hybrid=IH)
    f.CF = scalar
    f.CF.current_time = 3
    f.discretize()
    R = f.reconstruct(xi, et, sg, ravel=True)[1]
    RR = f.reconstruct(xi, et, sg, ravel=True, vectorized=True, value_only=True)[0]
    MR = f.do.make_reconstruction_matrix_on_grid(xi, et, sg)
    for _, i in enumerate(MR):
        r = MR[i] @ f.cochain.local[i]
        np.testing.assert_almost_equal(np.max(np.abs(r - R[i][0])), 0)
        np.testing.assert_almost_equal(np.max(np.abs(r - RR[_,:])), 0)

    LnEnF = f.do.compute_Ln_energy(vectorized=False)
    LnEnT = f.do.compute_Ln_energy(vectorized=True)
    np.testing.assert_almost_equal(LnEnF, LnEnT, decimal=5)

    #------------ 3-form -------------------------------------------------------
    f = FC('3-f', hybrid=IH)
    f.CF = scalar
    f.CF.current_time = 4
    f.discretize()
    R = f.reconstruct(xi, et, sg, ravel=True)[1]
    RR = f.reconstruct(xi, et, sg, ravel=True, vectorized=True, value_only=True)[0]
    MR = f.do.make_reconstruction_matrix_on_grid(xi, et, sg)
    for _, i in enumerate(MR):
        r = MR[i] @ f.cochain.local[i]
        np.testing.assert_almost_equal(np.max(np.abs(r - R[i][0])), 0)
        np.testing.assert_almost_equal(np.max(np.abs(r - RR[_,:])), 0)

    return 1


def test_Form_NO12_weak_curl():
    """"""
    if RANK == MASTER_RANK:
        print(f"~~~ [test_Form_NO12_weak_curl]... ", flush=True)
    else:
        pass

    # --- test 1: non-hybrid; periodic boundary condition ----------------------------------
    def u(t, x, y, z): return - np.sin(2*np.pi*x) * np.cos(2*np.pi*y) * np.sin(2*np.pi*z) + t/1.554
    def v(t, x, y, z): return np.cos(2*np.pi*x) * np.sin(2*np.pi*y) * np.cos(2*np.pi*z) + t*1.23
    def w(t, x, y, z): return np.sin(2*np.pi*x) * np.sin(2*np.pi*y) * np.cos(2*np.pi*z) + t/2.196
    mesh = MeshGenerator('crazy_periodic', c=0.0)([8, 10, 6], EDM=None)
    space = SpaceInvoker('polynomials')([3,2,4])
    FC = FormCaller(mesh, space)
    U = FC('vector', (u, v, w))
    curl_U = U.numerical.curl

    w1 = FC('1-f', hybrid=False)  # w1 = curl (u2)
    u2 = FC('2-f', hybrid=False)  # w1 = curl (u2)

    u2.CF = U
    u2.CF.current_time = 1
    u2.discretize()

    M1 = w1.matrices.mass

    M2 = u2.matrices.mass
    E12 = w1.matrices.incidence.T
    b = concatenate([E12 @ M2 @ u2.cochain.EWC,])
    A = bmat(([M1,],))

    A.gathering_matrices = (w1, w1)
    b.gathering_matrix = w1
    LS = LinearSystem(A, b)

    result = LS.solve('GMRES')(0, restart=50, maxiter=10)

    w1.CF = curl_U
    w1.CF.current_time = 1
    result[0].do.distributed_to(w1)

    assert w1.error.L() < 0.06, f"{w1.error.L()}! periodic boundary condition test fails."

    # --- test 2:  0 tangent velocity boundary condition ----------------------
    def u(t, x, y, z): return - np.cos(1.11*np.pi*x) * np.sin(2*np.pi*y) * np.sin(2*np.pi*z) + t/1.554
    def v(t, x, y, z): return np.sin(2*np.pi*x) * np.cos(0.85*np.pi*y) * np.sin(2*np.pi*z) + t*1.23
    def w(t, x, y, z): return np.sin(2*np.pi*x) * np.sin(2*np.pi*y) * np.cos(1.3*np.pi*z) + t/2.196
    mesh = MeshGenerator('crazy', c=0.0)([8, 10, 6], EDM=None)
    space = SpaceInvoker('polynomials')([3,2,4])
    FC = FormCaller(mesh, space)
    U = FC('vector', (u, v, w))
    curl_U = U.numerical.curl

    w1 = FC('1-f', hybrid=False)  # w1 = curl (u2)
    u2 = FC('2-f', hybrid=False)  # w1 = curl (u2)

    u2.CF = U
    u2.CF.current_time = 0
    u2.discretize()

    M1 = w1.matrices.mass

    M2 = u2.matrices.mass
    E12 = w1.matrices.incidence.T
    b = concatenate([E12 @ M2 @ u2.cochain.EWC,])
    A = bmat(([M1,],))

    A.gathering_matrices = (w1, w1)
    b.gathering_matrix = w1
    LS = LinearSystem(A, b)

    result = LS.solve('GMRES')(0, restart=50, maxiter=10)

    w1.CF = curl_U
    w1.CF.current_time = 0
    result[0].do.distributed_to(w1)

    assert w1.error.L() < 0.06, f"0 tangent velocity boundary condition test fails."

    return 1

if __name__ == '__main__':
    # mpiexec -n 4 python tests/objects/CSCG/_3d/unittests/standard_forms/general.py
    # test_Form_NO4_cross_product_1()
    # test_Form_NO5_cross_product_2()
    # test_Form_NOx_cross_product_3()
    # test_Form_NOx1_cross_product_4()
    test_Form_NO12_weak_curl()
    test_Form_NO2_mass_matrix()