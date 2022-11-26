# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/22 6:32 PM
"""
from collections.abc import Iterable


def single_list(list, ignore_types=(str)):
    for item in list:
        if isinstance(item, Iterable) and not isinstance(item, ignore_types):
            yield from single_list(item, ignore_types=(str))
        else:
            yield item