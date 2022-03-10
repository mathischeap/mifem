# mifem

v3.1.1

*mimetic finite element method*

*mifem* (mimetic finite element method) is an open-source Python implementation of the mimetic spectral element method
(MSEM) coded by Yi. It is a research-oriented project in the sense that it aims at providing a fine base for quick
implementations of new ideas for the MSEM rather than struggling for a higher computational efficiency and a wider range of 
applications.
Nevertheless, significant efforts have be made to boost the computational efficiency of *mifem*.
The most obvious feature as a result of such efforts is the usage of CPU parallelization powered by the Python binding of the
Message Passing Interface (MPI) standard, the *mpi4py* package.

Currently, no detailed documentation is applied as the number of users is quite low. 

Do
> mpiexec python TESTS\mpi_unittest.py

to run all the tests.

by Yi Zhang | www.mathischeap.com