# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 12/2/2022 3:49 PM
"""
import sys

if './' not in sys.path:
    sys.path.append('./')
from components.freeze.main import FrozenOnly
import numpy as np
from scipy.sparse import csr_matrix


def ___Pr_gathering_matrix_checker___(gathering, indices):
    """

    Parameters
    ----------
    gathering
    indices

    Returns
    -------
    num_dofs:

    """
    MAX = list()
    MIN = list()
    SET = set()
    for i in indices:
        gv = gathering[i]

        MAX.append(np.max(gv))
        MIN.append(np.min(gv))

        SET.update(gv)

    MAX = max(MAX)
    MIN = min(MIN)
    assert MIN == 0, f"Gathering matrix must number from 0."
    num_dofs = int(MAX + 1)
    assert len(SET) == num_dofs, f"gathering matrix wrong, it misses some numbers."

    return num_dofs


class VectorAssembler(FrozenOnly):
    """A local assembler. It does not assemble data across cores.

    It is costing, do not use it for large vectors.
    """

    def __init__(self, gathering):
        """

        Parameters
        ----------
        gathering :

        """
        if isinstance(gathering, dict):
            self._indices_ = gathering.keys()

        else:
            self._indices_ = range(len(gathering))

        self._gathering_ = gathering
        self._default_routine_ = 'basic'
        self._num_dofs_ = ___Pr_gathering_matrix_checker___(gathering, self._indices_)
        self._LEN_ = len(gathering)
        self._freeze_self_()

    def __call__(self, TBA, mode, routine=None, **kwargs):
        """The output will be a 1d array instead of sparse column vector for example."""
        if routine is None:
            routine = self._default_routine_
        else:
            pass
        return getattr(self, '___' + routine + "___")(TBA, mode, **kwargs)


    def ___basic___(self, TBA, mode, scheme=1, dtype=None):
        """"""
        if isinstance(self._gathering_, dict):
            assert isinstance(TBA, dict), f"data must be in a dict"
        else:
            pass

        assert len(TBA) == self._LEN_, f"data shape wrong."

        if scheme == 1:

            GM = self._gathering_
            indices = self._indices_

            if dtype is None:
                output = np.zeros(self._num_dofs_)
            else:
                output = np.zeros(self._num_dofs_, dtype=dtype)

            if mode == 'replace':
                for i in indices:
                    output[GM[i]] = TBA[i]
            elif mode == 'add':
                for i in indices:
                    output[GM[i]] += TBA[i]
            else:
                raise NotImplementedError(f"mode={mode} is not implemented "
                                          f"for scheme {scheme} of routine: basic")
            return output

        else:
            raise NotImplementedError(f"scheme = {scheme} is not implemented for routine: basic.")


class MatrixAssembler(FrozenOnly):
    """"""

    def __init__(self, G0, G1):
        """

        Parameters
        ----------
        G0
        G1
        """
        if isinstance(G0, dict):
            self._indices_ = G0.keys()
        else:
            self._indices_ = range(len(G0))

        self._G0_ = G0
        self._G1_ = G1
        self._default_routine_ = 'basic'
        self._num_dofs_0_ = ___Pr_gathering_matrix_checker___(G0, self._indices_)
        self._num_dofs_1_ = ___Pr_gathering_matrix_checker___(G1, self._indices_)
        self._LEN_ = len(G0)
        assert self._LEN_ == len(G1), f"gathering matrices length does not match."
        self._freeze_self_()

    def __call__(self, TBA, mode, routine=None, format='array', **kwargs):
        """

        Parameters
        ----------
        TBA :
            The data to be assembled.
        mode
        routine
        format
        kwargs

        Returns
        -------

        """
        if routine is None:
            routine = self._default_routine_
        else:
            pass
        return getattr(self, '___' + routine + "___")(TBA, mode, format=format, **kwargs)

    def ___basic___(self, TBA, mode, format='array', scheme=1):
        """"""
        if isinstance(self._G0_, dict):
            assert isinstance(TBA, dict), f"data must be in a dict"
        else:
            pass
        assert len(TBA) == self._LEN_, f"data shape wrong."

        if scheme == 1:

            output = np.zeros((self._num_dofs_0_, self._num_dofs_1_))

            for i in self._indices_:

                g0 = self._G0_[i]
                g1 = self._G1_[i]

                data = TBA[i]

                if mode == 'replace':
                    for j, m in enumerate(g0):
                        output[m, g1] = data[j, :]
                elif mode == 'add':
                    for j, m in enumerate(g0):
                        output[m, g1] += data[j, :]
                else:
                    raise NotImplementedError(
                        f"mode={mode} is not implemented for basic matrix assemble."
                    )

            if format == 'array':
                return output
            elif format == 'csc':
                return csr_matrix(output).T  # automatically a column csc matrix
            elif format == 'csr':
                return csr_matrix(output)   # a row csr matrix.
            else:
                raise NotImplementedError(f"cannot return {format} assembled matrix.")

        else:
            raise NotImplementedError(f"Not implemented for assembler scheme: {scheme}.")


if __name__ == '__main__':
    # mpiexec -n 1 python components/assemblers.py
    G = {1: [0, 1, 2],
         2: [1, 2, 3],
         3: [4, 5, 6, 0, 3, 1, 2],
         '4': np.array([3, 4, 7, 8, 9])}
    D = {
        1: [1, 1, 1],
        2: [1, 1, 1],
        3: [1, 1, 1, 1, 1, 1, 1],
        '4': [1, 1, 1, 1, 1]
    }

    output = VectorAssembler(G)(D, 'add')
    print(output)

    G = [[0, 1, 2],
         [1, 2, 3],
         [4, 5, 6, 0, 3, 1, 2],
         np.array([3, 4, 7, 8, 9])]
    D = (
        [1, 1, 1],
        [1, 1, 1],
        [1, 1, 1, 1, 1, 1, 1],
        [1, 1, 1, 1, 1]
    )

    output = VectorAssembler(G)(D, 'add')
    print(output)
