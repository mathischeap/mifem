# -*- coding: utf-8 -*-
"""
all unittests for 3d CSCG edge forms.
"""
import sys
if './' not in sys.path:
    sys.path.append('./')

from root.config.main import *
from root.save import save, read
import random
import os
from tests.objects.CSCG._3d.randObj.form_caller import random_FormCaller_of_total_load_around
from objects.CSCG._3d.master import ExactSolutionSelector, FormCaller, MeshGenerator, SpaceInvoker


def test_edge_forms_No0_save_read():
    """"""
    if RANK == MASTER_RANK:
        print("-ef- [test_edge_forms_No0_save_read] ...... ", flush=True)

    def p(t, x, y, z):
        return - 6 * np.pi * \
            np.sin(2.5*np.pi*x) * \
            np.sin(2.111*np.pi*y) * \
            np.sin(1.34*np.pi*z) \
            + 0.55 * t

    # # ------- below code makes the test file: `___unittest_edge_save_read_000___.mi`
    # mesh = MeshGenerator('crazy', c=0.0)([6, 5, 7], EDM='chaotic')
    # space = SpaceInvoker('polynomials')([('Lobatto', 2), ('Lobatto', 3), ('Lobatto', 1)])
    # FC = FormCaller(mesh, space)
    # scalar = FC('scalar', p)
    # e0 = FC('0-e')
    # e00 = FC('0-e')
    # e0.CF = scalar
    # e0.CF.current_time = 5.233
    # e0.discretize()
    #
    # ESS = ExactSolutionSelector(mesh)('icpsNS:sincosRD')
    # e00.CF = ESS.pressure
    # e00.CF.current_time = 2.875
    # e00.discretize()
    #
    # save([e0, e00], '___unittest_edge_save_read_000___.mi')
    # # ======================== do NOT DELETE ABOVE CODES ====================================

    absolute_path = os.path.dirname(__file__)
    E0, E00 = read(absolute_path + '/auxiliaries/___unittest_edge_save_read_000___.mi')

    if 100 in E0.cochain.local:
        np.testing.assert_array_almost_equal(
            E0.cochain.local[100],
            np.array([10.16985183,  0.97463385, -7.43387363, -0.12659664,  3.66254627,
                      7.12750344, 10.48126306,  0.89333921, -7.87427561, -0.25492221,
                      3.69604597,  7.30898322, 10.16985183,  7.92343458,  3.03542396,
                      -0.12659664, -7.43387363, -4.25695988,  2.65573103,  7.12750344,
                      10.48126306,  8.13890668,  3.04214076, -0.25492221, -7.87427561,
                      -4.56168344,  2.64623204,  7.30898322, 10.16985183, 10.48126306,
                      -7.43387363, -7.87427561, -0.12659664, -0.25492221,  7.12750344,
                      7.30898322])
        )
        np.testing.assert_array_almost_equal(
            E00.cochain.local[100],
            np.array(
                [
                    0.96194596,  0.69644976,  0.24434041,  0.03739119, -0.46726863,
                    -0.8467242,  0.38613605, -0.12681729, -0.60579004, -0.75797172,
                    -0.98256647, -0.94388333,  0.96194596,  0.81149797,  0.37530504,
                    0.03739119,  0.24434041, -0.10031749, -0.6150676, -0.8467242,
                    0.38613605,  0.04909339, -0.49068179, -0.75797172, -0.60579004,
                    -0.84043445, -0.99994247, -0.94388333,  0.96194596,  0.38613605,
                    0.24434041, -0.60579004,  0.03739119, -0.75797172, -0.8467242,
                    -0.94388333
                 ]
            )
        )

    if 125 in E0.cochain.local:
        np.testing.assert_array_almost_equal(
            E0.cochain.local[125],
            np.array(
                [
                    0.06165154, -1.07090409, -0.57134205,  7.3344841,  9.12644186,
                    8.33602233, -2.60671023, -4.81225353, -3.83940443, 11.55643261,
                    15.04609823, 13.50683212,  0.06165154,  2.14687699,  5.599493,
                    7.3344841, -0.57134205,  1.98252713,  6.21110088,  8.33602233,
                    -2.60671023,  1.45406606,  8.1777041, 11.55643261, -3.83940443,
                    1.1340105,  9.3687517, 13.50683212,  0.06165154, -2.60671023,
                    -0.57134205, -3.83940443,  7.3344841, 11.55643261,  8.33602233,
                    13.50683212
                ]
            )
        )
        np.testing.assert_array_almost_equal(
            E00.cochain.local[125],
            np.array(
                [
                    -0.738119, -0.97656508, -0.95334133, -0.86974113, -0.50646381,
                    -0.00747991, -0.98768834, -0.77714596, -0.35836795, -0.15643447,
                    0.35836795,  0.77714596, -0.738119, -0.92369066, -0.98578828,
                    -0.86974113, -0.95334133, -0.79365351, -0.34740838, -0.00747991,
                    -0.98768834, -0.87546191, -0.48328713, -0.15643447, -0.35836795,
                    -0.01919202,  0.51652869,  0.77714596, -0.738119, -0.98768834,
                    -0.95334133, -0.35836795, -0.86974113, -0.15643447, -0.00747991,
                    0.77714596
                ]
                     )
        )

    assert E0.mesh._EDM_ == 'chaotic' and E00.mesh._EDM_ == 'chaotic'

    # --- read time save and read -------------------------------------
    if RANK == MASTER_RANK:
        load = random.randint(50, 199)
        t = random.random()
    else:
        load, t = None, None
    load, t = COMM.bcast([load, t], root=MASTER_RANK)
    FC = random_FormCaller_of_total_load_around(load, exclude_periodic=True)
    ESS = ExactSolutionSelector(FC.mesh)('icpsNS:sincosRD')
    scalar = FC('scalar', p)

    e0 = FC('0-e', orientation='outer', numbering_parameters={'scheme_name': 'Naive'}, name='test_0edge')
    e0.CF = scalar
    e0.CF.current_time = t
    e0.discretize()
    save(e0, '___unittest_edge_e0_save___.mi')
    E0 = read('___unittest_edge_e0_save___.mi')
    c0 = e0.cochain.local
    C0 = E0.cochain.local
    for mei in c0:
        np.testing.assert_array_almost_equal(c0[mei], C0[mei])
    assert E0.standard_properties.name == 'test_0edge'
    assert E0.numbering._scheme_name_ == 'Naive'
    assert E0.numbering._parameters_ == dict()
    assert E0.orientation == 'outer'

    e0 = FC('0-e', orientation='inner', numbering_parameters={'scheme_name': 'Naive'}, name='test_0edge111')
    e0.CF = ESS.pressure
    e0.CF.current_time = t + 10
    e0.discretize()
    save(e0, '___unittest_edge_e0_save1___.mi')
    E0 = read('___unittest_edge_e0_save1___.mi')
    c0 = e0.cochain.local
    C0 = E0.cochain.local
    for mei in c0:
        np.testing.assert_array_almost_equal(c0[mei], C0[mei])
    assert E0.standard_properties.name == 'test_0edge111'
    assert E0.numbering._scheme_name_ == 'Naive'
    assert E0.numbering._parameters_ == dict()
    assert E0.orientation == 'inner'

    if RANK == MASTER_RANK:
        os.remove('___unittest_edge_e0_save___.mi')
        os.remove('___unittest_edge_e0_save1___.mi')

    return 1


def test_edge_forms_No1_0edge_Rd_and_Rc():
    """"""
    if RANK == MASTER_RANK:
        print("-ef- [test_edge_forms_No1_0edge_Rd_and_Rc] ...... ", flush=True)

    mesh = MeshGenerator('crazy', c=0.095)([3, 4, 5])
    space = SpaceInvoker('polynomials')([('Lobatto', 5), ('Lobatto', 4), ('Lobatto', 3)])
    FC = FormCaller(mesh, space)
    e0 = FC('0-e')
    def p(t, x, y, z): return - 6 * np.pi * np.sin(2*np.pi*x) * np.sin(2*np.pi*y) * np.sin(2*np.pi*z) + 1.2 * t
    scalar = FC('scalar', p)
    e0.CF = scalar
    e0.CF.current_time = 0
    e0.discretize()
    assert e0.error.L(n='infinity', quad_density=5000000) < 0.1

    return 1


if __name__ == '__main__':
    # mpiexec -n 6 python objects\CSCG\_3d\__tests__\unittests\edge_forms.py
    test_edge_forms_No0_save_read()
