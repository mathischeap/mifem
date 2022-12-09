# -*- coding: utf-8 -*-

def break_list_into_parts(LIST, parts):
    """For example:

    .. doctest::

        >>> break_list_into_parts([1,2,3,4,5], [3,2])
        ([1, 2, 3], [4, 5])

    :param LIST:
    :param parts:
    :return:
    """
    start = 0
    RET = tuple()
    for p in parts:
        RET += (LIST[start:start+p],)
        start += p
    return RET