# -*- coding: utf-8 -*-
from components.freeze.base import FrozenOnly

class SchurSolverBase(FrozenOnly):
    """A base for all Schur solvers."""
    def __init__(self, rank, blocks, routine, name):
        """

        :param routine: Which particular routine we are going to use?
        :param name: The name of this solving process.
        """
        self._rank_ = rank
        self._blocks_ = blocks
        self._routine_ = routine
        self._name_ = name
        self._freeze_self_()


    def __call__(self, *args, **kwargs):
        """
        :return: Return a tuple of 5 outputs:

            1. (LocallyFullVector) results -- The result vector.
            2. (int) info -- The info which provides convergence information:

                * 0 : successful exit
                * >0 : convergence to tolerance not achieved, number of iterations
                * -1 : divergence

            3. (float) beta -- The residual.
            4. (int) ITER -- The number of outer iterations.
            5. (str) message

        """
        raise NotImplementedError()