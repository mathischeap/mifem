# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/23 9:20 PM
"""
import sys

if './' not in sys.path: sys.path.append('/')

from objects.mpRfT._2d.__tests__.Random.mesh import rm
from objects.mpRfT._2d.master import FormCaller


def rf(*args, **kwargs):
    mesh = rm(*args, **kwargs)
    return FormCaller(mesh)


if __name__ == '__main__':
    # mpiexec -n 4 python objects/mpRfT/_2d/__tests__/Random/form.py
    fc = rf(100)

    f = fc('0-f-o')

    print(f)