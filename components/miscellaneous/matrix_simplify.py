# -*- coding: utf-8 -*-
import numpy as np


def rsmat(arbmat):
    """ Convert an arbitrary matrix to its simplest format """
    arbmat = arbmat.astype(float)
    row_number, column_number = arbmat.shape
    if row_number == 1:
        if arbmat[0, 0] != 0:
            return arbmat / arbmat[0, 0]
        else:
            return arbmat
    else:
        rc_number = min(row_number, column_number)
        anarbmat = arbmat.copy()
        r = 0
        for n in range(rc_number):
            s_row = -1
            for i in arbmat[r:row_number, n]:
                s_row += 1
                if abs(i) > 1e-10:
                    anarbmat[r, :] = arbmat[s_row + r, :]
                    for j in range(r, row_number):
                        if j < s_row + r:
                            anarbmat[j + 1, :] = arbmat[j, :]
                    arbmat = anarbmat.copy()
            if abs(anarbmat[r, n]) > 1e-10:
                anarbmat[r, :] = anarbmat[r, :] / anarbmat[r, n]
                for i in range(row_number):
                    if i != r:
                        anarbmat[i, :] -= \
                            anarbmat[i, n] * anarbmat[r, :]
            arbmat = anarbmat.copy()
            if abs(arbmat[r, n]) < 1e-10:
                r = r
            else:
                r = r + 1
        for m in range(column_number):
            if abs(arbmat[-1, m]) > 1e-10:
                arbmat[-1, :] = arbmat[-1, :] / arbmat[-1, m]
                for i in range(row_number - 1):
                    arbmat[i, :] -= \
                        arbmat[i, m] * arbmat[-1, :]
                break

        return arbmat


if __name__ == '__main__':

    a = np.array([
            [1, 0, 2, -1, 1, 3],
            [2, 0, 3, 1, 4, 5],
            [2, 0, 3, 1, 4, 5]])


    a1 = rsmat(a)

    print(a1)