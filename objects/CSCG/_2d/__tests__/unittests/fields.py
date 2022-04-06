
import sys
if './' not in sys.path: sys.path.append('./')
from root.config.main import *
import random
from objects.CSCG._2d.__tests__.random_objects.form_caller import random_FormCaller_of_total_load_around

def test_Fields_NO1_vector():
    """"""
    if rAnk == mAster_rank:
        load = random.randint(100,1000)
        print(f"~~~ [test_Fields_NO1_vector] @ load= {load}... ", flush=True)
    else:
        load = None

    load = cOmm.bcast(load, root=mAster_rank)
    FC = random_FormCaller_of_total_load_around(load)

    I, J = random.randint(2,10), random.randint(4,8)
    x = np.linspace(-0.9-random.random()/10, 0.9+random.random()/10, I)
    y = np.linspace(-0.8-random.random()/5, 0.8+random.random()/5, J)
    x, y = np.meshgrid(x, y, indexing='ij')
    t = random.random()

    def u(t, x, y): return - np.pi * np.sin(2.56*np.pi*x) * np.cos(3.12*np.pi*y) * np.sin(t)/1.554
    def v(t, x, y): return np.pi * np.cos(2.112*np.pi*x) * np.sin(1.98*np.pi*y) * np.sin(t)*1.23
    def du_dt(t, x, y): return - np.pi * np.sin(2.56*np.pi*x) * np.cos(3.12*np.pi*y) * np.cos(t)/1.554
    def dv_dt(t, x, y): return np.pi * np.cos(2.112*np.pi*x) * np.sin(1.98*np.pi*y) * np.cos(t)*1.23

    v = FC('vector', (u, v))
    # test time derivatives ---------------------------------------------------
    v_t = v.numerical.time_derivative
    v_t.current_time = t
    R_xyz, R_v = v_t.reconstruct(x, y)
    for i in R_xyz:
        assert np.max(np.abs(R_v[i][0] - du_dt(t, *R_xyz[i]))) < 1e-7
        assert np.max(np.abs(R_v[i][1] - dv_dt(t, *R_xyz[i]))) < 1e-7

    #  ---- test numerical curl -------------------------------------------------------------
    def dv_dx__m__du_dy(t, x, y):
        return -2.112*np.pi**2 * np.sin(2.112*np.pi*x) * np.sin(1.98*np.pi*y) * np.sin(t)*1.23 - \
               3.12 * np.pi**2 * np.sin(2.56 * np.pi * x) * np.sin(3.12 * np.pi * y) * np.sin(t)/1.554
    curl_v = v.numerical.curl
    curl_v.current_time = t
    R_xyz, R_v = curl_v.reconstruct(x, y)
    for i in R_xyz:
        assert np.max(np.abs(R_v[i][0] - dv_dx__m__du_dy(t, *R_xyz[i]))) < 1e-7


    #  ---- test numerical div -------------------------------------------------------------
    def du_dx__p__dv_dy(t, x, y):
        return - 2.56*np.pi**2 * np.cos(2.56*np.pi*x) * np.cos(3.12*np.pi*y) * np.sin(t)/1.554 + \
               1.98*np.pi**2 * np.cos(2.112*np.pi*x) * np.cos(1.98*np.pi*y) * np.sin(t)*1.23
    div_v = v.numerical.div
    div_v.current_time = t
    R_xyz, R_v = div_v.reconstruct(x, y)
    for i in R_xyz:
        assert np.max(np.abs(R_v[i][0] - du_dx__p__dv_dy(t, *R_xyz[i]))) < 1e-7

    return 1


def test_Fields_NO2_scalar():
    """"""
    if rAnk == mAster_rank:
        load = random.randint(100,1000)
        print(f"~~~ [test_Fields_NO2_scalar] @ load= {load}... ", flush=True)
    else:
        load = None

    load = cOmm.bcast(load, root=mAster_rank)
    FC = random_FormCaller_of_total_load_around(load)

    I, J = random.randint(2,10), random.randint(4,8)
    x = np.linspace(-0.9-random.random()/10, 0.9+random.random()/10, I)
    y = np.linspace(-0.8-random.random()/5, 0.8+random.random()/5, J)
    x, y = np.meshgrid(x, y, indexing='ij')
    t = random.random()
    def f(t, x, y): return - np.pi * np.sin(2.56*np.pi*x) * np.cos(3.12*np.pi*y) * np.sin(t) / 1.554
    def df_dx(t, x, y): return - 2.56* np.pi**2 * np.cos(2.56*np.pi*x) * np.cos(3.12*np.pi*y) * np.sin(t) / 1.554
    def df_dy(t, x, y): return 3.12* np.pi**2 * np.sin(2.56*np.pi*x) * np.sin(3.12*np.pi*y) * np.sin(t) / 1.554
    def df_dt(t, x, y): return - np.pi * np.sin(2.56*np.pi*x) * np.cos(3.12*np.pi*y) * np.cos(t) / 1.554
    w = FC('scalar', f)

    # test time derivatives ---------------------------------------------------
    w_t = w.numerical.time_derivative
    w_t.current_time = t
    R_xyz, R_v = w_t.reconstruct(x, y)
    for i in R_xyz:
        assert np.max(np.abs(R_v[i][0] - df_dt(t, *R_xyz[i]))) < 1e-7


    #----------- test numerical gradient ------------------------------------------
    grad_w = w.numerical.grad
    grad_w.current_time = t
    R_xyz, R_v = grad_w.reconstruct(x, y)
    for i in R_xyz:
        assert np.max(np.abs(R_v[i][0] - df_dx(t, *R_xyz[i]))) < 1e-7
        assert np.max(np.abs(R_v[i][1] - df_dy(t, *R_xyz[i]))) < 1e-7

    #----------- test numerical curl ------------------------------------------
    curl_w = w.numerical.curl
    curl_w.current_time = t
    R_xyz, R_v = curl_w.reconstruct(x, y)
    for i in R_xyz:
        assert np.max(np.abs(R_v[i][0] - df_dy(t, *R_xyz[i]))) < 1e-7
        assert np.max(np.abs(R_v[i][1] + df_dx(t, *R_xyz[i]))) < 1e-7


    return 1








if __name__ == '__main__':
    # mpiexec -n 4 python objects\CSCG\_2d\__tests__\unittests\fields.py
    test_Fields_NO1_vector()
    test_Fields_NO2_scalar()