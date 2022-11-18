# mifem

v3.2.2

*mimetic finite element method*

*mifem* is an open-source Python implementation of the mimetic spectral element method
(MSEM) coded by Yi. It is a research-oriented project in the sense that it aims at providing a fine base for quick
implementations of new ideas for the MSEM rather than struggling for a higher computational efficiency and a wider range of 
applications.
Nevertheless, significant efforts have be made to boost the computational efficiency of *mifem*.
The most obvious feature as a result of such efforts is the usage of CPU parallelization powered by the Python binding of the
Message Passing Interface (MPI) standard, the *mpi4py* package.


To start up, please first check the installation of the dependent packages in *requirements.txt*. 
Then you can direct to the dir where mifem library locates and run, for example,
```
$ mpiexec python tests/run.py
```
This will do all the (>100) tests to fully validate the code in your machine.


by Yi Zhang | www.mathischeap.com