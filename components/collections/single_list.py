# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/22 6:32 PM
"""
from collections.abc import Iterable


def single_list(List):
    """Iterate over all levels except string.

    Examples
    --------
        >>> List = [[1, (2, 2.25, 2.75), 3, [3.5, 4], 'abc'], 5, 6]
        >>> for i in single_list(List): print(i)
        1
        2
        2.25
        2.75
        3
        3.5
        4
        abc
        5
        6

    Parameters
    ----------
    List

    Returns
    -------

    """
    for item in List:
        if isinstance(item, Iterable) and not isinstance(item, str):
            yield from single_list(item)
        else:
            yield item


if __name__ == '__main__':

    from doctest import testmod

    testmod()
