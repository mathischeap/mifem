# -*- coding: utf-8 -*-


def ___split_A___(indptr, splitting_factor, PS):
    """A private function to help split global matrix into parts.

    :param indptr:
    :param splitting_factor:
    :param PS:
    :return:
    """
    assert len(indptr) == PS+1
    ST = 0
    __ = splitting_factor
    SL = [0, ]
    for i, ind in enumerate(indptr):
        if ind >= __:
            ST += 1
            SL.append(i)
            __ += splitting_factor
    if SL[-1] != PS:
        ST += 1
        SL.append(PS)
    assert ST >= 1 and ST == len(SL)-1
    return ST, SL
