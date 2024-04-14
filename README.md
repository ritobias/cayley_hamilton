# Iterative Cayley-Hamilton
header-only C++ implementation of the iterative Cayley-Hamilton method for the evaluation of matrix power series (and their differentials) of complex square matrices of arbitrary size. For a more detailed description of the method, please see [arXiv:2404.07704](https://arxiv.org/abs/2404.07704).

## Content

### cayley_hamilton.h
contains the implementation of a callable cayley_hamilton class with three overloads, computing
- only the matrix polynomial,
- the matrix polynomial and its derivative in a given direction,
- the matrix polynomial and all its derivatives.

See comments in the header file for more details.


### sun_algebra.h
contains the implementation of a sun_algebra class, which provides:
- hermitian su(n) basis as sparse_mat objects
- member functions to transform between su(n) vectors and su(n) matrices in either anti-hermition or hermitian rep.:
   * anti-hermitian: get_alg_mat_ah(invec,outmat) <--> proj_ah(inmat,outvec)
   * hermitiaon: get_alg_mat_h(invec,outmat) <--> proj_h(inmat,outvec)
- member function to take matrix log of SU(N) matrices and return corresponding su(n) vector: log_ah(inmat,outvec)
- member function to compute the SU(N) matrix corresponding to an su(n) vector: get_grp_mat(invec,outmat) (using the cayley_hamilton class)

### test.cpp
CPP file to run some basic tests with the cayley_hamilton and sun_algebra classes.
