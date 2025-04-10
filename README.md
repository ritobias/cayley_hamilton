# Iterative Cayley-Hamilton
header-only C++ implementation of the iterative Cayley-Hamilton method for the evaluation of matrix power series (and their differentials) of complex square matrices of arbitrary size. For a more detailed description of the method, please see [arXiv:2404.07704](https://arxiv.org/abs/2404.07704).

## Content

### cayley_hamilton.h
#### cayley_hamilton class
implementation of the iterative Cayley-Hamilton method for computing matrix power series. 
Callable class has three overloads, computing
- only the matrix polynomial,
- the matrix polynomial and its derivative in a given direction,
- the matrix polynomial and all its derivatives.

has additional subroutines:
- get_r_k computes and returns the Cayley-Hamilton decomposition coefficents for the matrix polynomial and of its derivatives
- get_r_dr computes and returns the Cayley-Hamilton decomposition coefficients of the matrix polynomial and the derivatives of these coefficients 
- ch_mult computes the Cayley-Hamilton decomposition coefficients of the product of two matrix polynomials with given CH decomposition coefficients_

#### chexp class
implementation of the iterative Cayley-Hamilton method with scaling and squaring for computing matrix exponentials. 
Callable class has three overloads, computing
- only the matrix exponential,
- the matrix exponential and its derivative in a given direction,
- the matrix exponential and all its derivatives.

has additional subroutines:
- get_r_k computes and returns the Cayley-Hamilton decomposition coefficents for the matrix exponential and of its derivatives
- get_r_dr computes and returns the Cayley-Hamilton decomposition coefficients of the matrix exponential and the derivatives of these coefficients 


#### nvexp class
callable class, implementing the matrix exponentiation using naive Taylor series definition in combination with scaling and squaring. 


See comments in the header file for more details.


### sun_algebra.h
#### sun_algebra class
implementation of a sun_algebra class, which provides:
- hermitian su(n) basis as sparse_mat objects
- member functions to transform between su(n) vectors and su(n) matrices in either anti-hermition or hermitian rep.:
   * anti-hermitian: get_alg_mat_ah(invec,outmat) <--> proj_ah(inmat,outvec)
   * hermitiaon: get_alg_mat_h(invec,outmat) <--> proj_h(inmat,outvec)
- member function to take matrix log of SU(N) matrices and return corresponding su(n) vector: log_ah(inmat,outvec)
- member function to compute the SU(N) matrix corresponding to an su(n) vector: get_grp_mat(invec,outmat) (using the cayley_hamilton class)

### test.cpp
CPP file to run some basic tests with the cayley_hamilton, chexp, and sun_algebra classes.
