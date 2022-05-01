
import sys
if './' not in sys.path: sys.path.append('./')
from root.config.main import *

import random, time
import os

from screws.functions._3d_space.scaling import ScalingFunc
from screws.functions._3d_space.opposite import Opposite
from screws.functions._3d_space.constant import CFG
from screws.functions._3d_space._0_ import _0_
from screws.functions._3d_space.Cartesian_spherical_coordinate_switcher import CartSphSwitcher
from screws.functions._3d_space.Cartesian_cylinder_coordinate_switcher import CartCylSwitcher
from screws.numerical.time_plus_3d_space.partial_derivative_as_functions import NumericalPartialDerivative_txyz_Functions
from screws.numerical.time_plus_3d_space.partial_derivative import NumericalPartialDerivative_txyz
from screws.emails.plain import SendAdminAnEmail, SendAdminAnHTMLEmail
from screws.miscellaneous.generalized_piecewise_function import genpiecewise

from functools import partial

def SIN(x, y, z): return np.sin(x * y * z)


def test_SCREWS_NO1_3d_functions():
    """ """
    if rAnk == mAster_rank:
        print("--- [test_SCREWS_NO1_3d_functions] ...... ", flush=True)

    i = random.randint(2,10)
    j = random.randint(2,10)

    x = np.random.rand(i, j)
    y = np.random.rand(i, j)
    z = np.random.rand(i, j)
    A = SIN(x, y, z)
    f = random.uniform(1, 10)
    B = ScalingFunc(f)(SIN)(x, y, z)
    np.testing.assert_array_almost_equal(A*f, B)
    C = Opposite(SIN)()(x, y, z)
    np.testing.assert_array_almost_equal(-A, C)

    D = CFG(f)()(x, y, z)
    np.testing.assert_array_almost_equal(D, np.ones((i, j))*f)

    E = _0_(x, y, z)
    np.testing.assert_array_almost_equal(E, np.zeros((i, j)))

    a, b, c = CartCylSwitcher.cart2cyl(x, y, z)
    X, Y, Z = CartCylSwitcher.cyl2cart(a, b, c)
    np.testing.assert_array_almost_equal(x, X)
    np.testing.assert_array_almost_equal(y, Y)
    np.testing.assert_array_almost_equal(z, Z)

    a, b, c = CartSphSwitcher.cart2sph(x, y, z)
    X, Y, Z = CartSphSwitcher.sph2cart(a, b, c)
    np.testing.assert_array_almost_equal(x, X)
    np.testing.assert_array_almost_equal(y, Y)
    np.testing.assert_array_almost_equal(z, Z)

    return 1



def test_SCREWS_NO2_sending_an_email_to_admin(force_to_do=False):
    """ """
    if rAnk == mAster_rank:
        print("-@- [test_SCREWS_NO2_sending_an_email_to_admin] ...... ", flush=True)

        # noinspection PyBroadException
        try:
            ct = time.ctime()
            ct_format = time.strptime(ct, '%a %b %d %H:%M:%S %Y')
            ct_str = time.strftime("%Y %m %d %H %M %S", ct_format)
            Y, M, D, H, _, _ = ct_str.split(' ')
            Y, M, D, H, = int(Y), int(M), int(D), int(H)
            Current_day = (Y - 1) * 365 * 24 + (M - 1) * 30 * 24 + D * 24 + H

            absolute_path = os.path.dirname(__file__)

            with open(absolute_path + '/auxiliaries/___private_developer_email_admin_1___.txt', 'r') as f:
                plain_last_time_day = f.readline()
                time_format = time.strptime(plain_last_time_day, '%a %b %d %H:%M:%S %Y')
                time_str = time.strftime("%Y %m %d %H %M %S", time_format)
                Y, M, D, H, _, _ = time_str.split(' ')
                Y, M, D, H = int(Y), int(M), int(D), int(H)
                plain_last_time_day = (Y - 1) * 365 * 24 + (M - 1) * 30 * 24 + D * 24 + H

            with open(absolute_path + '/auxiliaries/___private_developer_email_admin_2___.txt', 'r') as f:
                email_last_time_day = f.readline()
                time_format = time.strptime(email_last_time_day, '%a %b %d %H:%M:%S %Y')
                time_str = time.strftime("%Y %m %d %H %M %S", time_format)
                Y, M, D, H, _, _ = time_str.split(' ')
                Y, M, D, H = int(Y), int(M), int(D), int(H)
                email_last_time_day = (Y - 1) * 365 * 24 + (M - 1) * 30 * 24 + D * 24 + H

            D1 = Current_day-plain_last_time_day
            D2 = Current_day-email_last_time_day

            i = random.random()
            if force_to_do or i > 0.98 or D1 > 24 * 14: # 2% chance to do this test or did not test for more than 14 days
                text = 'A plain unittest message. You see this when a successful unittest just passes.'
                code = SendAdminAnEmail()(text)
                if code == 1:
                    with open(absolute_path + '/auxiliaries/___private_developer_email_admin_1___.txt', 'w') as f:
                        f.write(time.ctime())
                    if force_to_do:
                        print("   ~ Plain test email sent. (forced to do)", flush=True)
                    elif i > 0.98:
                        print("   ~ Plain test email sent. (2% chance)", flush=True)
                    else:
                        print("   ~ Plain test email sent. (It's time!)", flush=True)
                elif code == 0:
                    print("   ~ Plain test email sent failed. But it is OK, it is not an error.")
                else:
                    raise Exception(f"unknown plain email sending code {code}.")
            else: # we do not do the test.
                print("   ~ Plain email test skipped.", flush=True)

            i = random.random()
            if force_to_do or i > 0.98 or D2 > 24 * 7: # 2% chance to do this test or did not test for more than 7 days
                text = 'A HTML unittest message. You see this when a successful unittest just passes.'
                code = SendAdminAnHTMLEmail()(text)
                if code == 1:
                    with open(absolute_path + '/auxiliaries/___private_developer_email_admin_2___.txt', 'w') as f:
                        f.write(time.ctime())
                    if force_to_do:
                        print("   ~ HTML test email sent. (forced to do)", flush=True)
                    elif i > 0.98:
                        print("   ~ HTML test email sent. (2% chance)", flush=True)
                    else:
                        print("   ~ HTML test email sent. (It's time!)", flush=True)
                elif code == 0:
                    print("   ~ HTML test email sent failed. But it is OK, it is not an error.")
                else:
                    raise Exception(f"unknown HTML email sending code {code}.")
            else: # we do not do the test.
                print("   ~ HTML email test skipped.", flush=True)

        except:

            pass

    return 1



def test_SCREWS_NO3_4d_functions():
    """ """
    if rAnk == mAster_rank:
        print("-4- [test_SCREWS_NO3_4d_functions] ...... ", flush=True)

    def func(t, x, y, z): return np.sin(np.pi*x) * np.sin(np.pi*y) * np.sin(np.pi*z) * t
    def APt(t, x, y, z): return np.sin(np.pi*x) * np.sin(np.pi*y) * np.sin(np.pi*z) + 0*t
    def APx(t, x, y, z): return np.pi*np.cos(np.pi*x) * np.sin(np.pi*y) * np.sin(np.pi*z) * t
    def APy(t, x, y, z): return np.pi*np.sin(np.pi*x) * np.cos(np.pi*y) * np.sin(np.pi*z) * t
    def APz(t, x, y, z): return np.pi*np.sin(np.pi*x) * np.sin(np.pi*y) * np.cos(np.pi*z) * t

    t = random.random() * 100
    I, J, K = random.randint(10,20), random.randint(10,20), random.randint(10,20)
    x = np.random.rand(I, J, K)
    y = np.random.rand(I, J, K)
    z = np.random.rand(I, J, K)
    NP = NumericalPartialDerivative_txyz(func, t, x, y, z)
    assert all(NP.check_total(APt, APx, APy, APz))
    NPD4F = NumericalPartialDerivative_txyz_Functions(func)
    NPt = NPD4F('t')
    NPx = NPD4F('x')
    NPy = NPD4F('y')
    NPz = NPD4F('z')
    assert np.sum(np.abs(APt(t,x,y,z) - NPt(t, x, y, z))) < 1e-4, f"Partial t error too big!"
    assert np.sum(np.abs(APx(t,x,y,z) - NPx(t, x, y, z))) < 1e-4, f"Partial x error too big!"
    assert np.sum(np.abs(APy(t,x,y,z) - NPy(t, x, y, z))) < 1e-4, f"Partial y error too big!"
    assert np.sum(np.abs(APz(t,x,y,z) - NPz(t, x, y, z))) < 1e-4, f"Partial z error too big!"

    NPt_t = partial(NPt, t)
    NPx_t = partial(NPx, t)
    NPy_t = partial(NPy, t)
    NPz_t = partial(NPz, t)
    assert np.sum(np.abs(APt(t,x,y,z) - NPt_t(x, y, z))) < 1e-4, f"Partial t error too big!"
    assert np.sum(np.abs(APx(t,x,y,z) - NPx_t(x, y, z))) < 1e-4, f"Partial x error too big!"
    assert np.sum(np.abs(APy(t,x,y,z) - NPy_t(x, y, z))) < 1e-4, f"Partial y error too big!"
    assert np.sum(np.abs(APz(t,x,y,z) - NPz_t(x, y, z))) < 1e-4, f"Partial z error too big!"

    assert np.sum(np.abs(NPt(t,x,y,z) - NPt_t(x, y, z))) < 1e-10, f"Partial t error too big!"
    assert np.sum(np.abs(NPx(t,x,y,z) - NPx_t(x, y, z))) < 1e-10, f"Partial x error too big!"
    assert np.sum(np.abs(NPy(t,x,y,z) - NPy_t(x, y, z))) < 1e-10, f"Partial y error too big!"
    assert np.sum(np.abs(NPz(t,x,y,z) - NPz_t(x, y, z))) < 1e-10, f"Partial z error too big!"



def test_SCREWS_NO4_generalized_piecewise_function():
    """ """
    if rAnk == mAster_rank:
        print("-4- [test_SCREWS_NO4_generalized_piecewise_function] ...... ", flush=True)

    def FUNC1(t, x, y, z): return - 1 - x * y - z * t
    def FUNC2(t, x, y, z): return + 2 + x + y * z - t
    def FUNC3(t, x, y, z): return - 3 - x + y * z - t
    def FUNC4(t, x, y, z): return + 4 - x * y + z + t

    x = np.random.rand(10,10)
    y = np.random.rand(10,10)
    z = np.random.rand(10,10)
    t = 1

    cond1, cond2 = y <= 0.5, (y > 0.5) & (y<=1)
    RESULT = genpiecewise([t, x, y, z], [cond1, cond2], [FUNC3, FUNC4])
    #
    assert np.all(FUNC3(t, x[cond1], y[cond1], z[cond1]) == RESULT[cond1])
    assert np.all(FUNC4(t, x[cond2], y[cond2], z[cond2]) == RESULT[cond2])


    cond1, cond2, cond3, cond4 = y <= 0.25, \
                                 (0.25 < y) & (y <= 0.5), \
                                 (0.5  < y) & (y <= 0.75), \
                                 0.75 < y
    RESULT = genpiecewise([t, x, y, z], [cond1, cond2, cond3, cond4], [FUNC1, FUNC2, FUNC3, FUNC4])
    assert np.all(FUNC1(t, x[cond1], y[cond1], z[cond1]) == RESULT[cond1])
    assert np.all(FUNC2(t, x[cond2], y[cond2], z[cond2]) == RESULT[cond2])
    assert np.all(FUNC3(t, x[cond3], y[cond3], z[cond3]) == RESULT[cond3])
    assert np.all(FUNC4(t, x[cond4], y[cond4], z[cond4]) == RESULT[cond4])



if __name__ == '__main__':
    # mpiexec python __tests__\unittests\screws_.py

    test_SCREWS_NO2_sending_an_email_to_admin()