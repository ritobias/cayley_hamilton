/*
	header-only C++ implementation of an su(n) Lie-algebra class
	Copyright (C) 2024  Tobias Rindlisbacher

	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program.  If not, see <https://www.gnu.org/licenses/>.

	Email: ritobias@gmx.ch
*/
#pragma once
#include "cayley_hamilton.h"

// note: type definitions (ftype, ctype, etc.) are done in cayley_hamilton.h 

template<class T>
class sa_entry {
	// structure to hold row index (ind1), column index (ind2), and corresponding value (val) of a matrix component  
public:
	sa_entry():ind1(0),ind2(0),val(0) {

	}
	sa_entry(int ti1,int ti2,T tval): ind1(ti1),ind2(ti2),val(tval) {

	}
	void set(int ti1,int ti2,T tval) {
		ind1=ti1;
		ind2=ti2;
		val=tval;
	}
	void clear() {
		ind1=0;
		ind2=0;
		val=0;
	}
	~sa_entry() {

	}
	int ind1;
	int ind2;
	T val;
};

template<class T>
class sparse_mat {
	// sparse matrix class which stores the non-zero matrix elements as list of sa_entry objects
public:
	sparse_mat():n(0),nelem(0),elem(0) {

	}
	sparse_mat(int tn): n(tn) {
		elem=new sa_entry<T>[n];
		nelem=0;
	}
	void init(int tn) {
		if(elem!=0) {
			delete[] elem;
			elem=0;
		}

		if(tn>0) {
			n=tn;
			elem=new sa_entry<T>[n];
		}
		nelem=0;
	}
	void push(int ti1,int ti2,T tval) {
		elem[nelem].set(ti1,ti2,tval);
		++nelem;
	}
	void pop() {
		if(nelem>0) {
			elem[nelem-1].clear();
			--nelem;
		}
	}
	~sparse_mat() {
		delete[] elem;
		nelem=0;
	}
	void print() {
		T** tmat=new T*[n];
		for(int i=0; i<n; ++i) {
			tmat[i]=new T[n]();
		}
		for(int i=0; i<nelem; ++i) {
			tmat[elem[i].ind1][elem[i].ind2]=elem[i].val;
		}
		std::cout<<std::endl;
		for(int i=0; i<n; ++i) {
			for(int j=0; j<n; ++j) {
				std::cout<<tmat[i][j]<<" ";
			}
			std::cout<<std::endl;
		}
		std::cout<<std::endl;

		for(int i=0; i<n; ++i) {
			delete[] tmat[i];
		}
		delete[] tmat;
	}
	int nelem;
	sa_entry<T>* elem;
private:
	int n;
};

class sun_algebra {
	// class providing:
	// - hermitian su(n) basis as sparse_mat objects
	// - member functions to transform between su(n) vectors and su(n) matrices in either anti-hermition or hermitian rep.:
	//    -- anti-hermitian: get_alg_mat_ah(invec[],outmat[][]) <--> proj_ah(inmat[][],outvec[])
	//    -- hermitiaon: get_alg_mat_h(invec[],outmat[][]) <--> proj_h(inmat[][],outvec[])
	// - member function to take matrix log of SU(N) matrices and return corresponding su(n) vector: log_ah(inmat[][],outvec[])
	// - member function to compute the SU(N) matrix corresponding to an su(n) vector: get_grp_mat(invec[],outmat[][])
public:
	sun_algebra(int tn): ch(tn),tgen(0),telem0(0),telem(0) {
		n=tn;
		ngen=n*n-1;
		generator=new sparse_mat<ctype>[ngen];
		int i,j;
		for(i=0; i<ngen; ++i) {
			generator[i].init(n);
		}
		// generate the hermitian su(n) algebra basis matrices \lambda[i][][] , i=0,...,n^2-1:
		// (lambda[i][][] is stored in sparse form in generator[i])

		// symmetric generators
		int igen=0;
		for(i=0; i<n-1; ++i) {
			for(j=i+1; j<n; ++j) {
				sparse_mat<ctype>& tgen=generator[igen];
				tgen.push(i,j,1);
				tgen.push(j,i,1);
				++igen;
			}
		}
		// anti-symmetric generators
		for(i=0; i<n-1; ++i) {
			for(j=i+1; j<n; ++j) {
				sparse_mat<ctype>& tgen=generator[igen];
				tgen.push(i,j,ctype(0,-1));
				tgen.push(j,i,ctype(0,1));
				++igen;
			}
		}
		// diagonal generators
		for(i=0; i<n-1; ++i) {
			sparse_mat<ctype>& tgen=generator[igen];
			for(j=0; j<=i; ++j) {
				tgen.push(j,j,std::sqrt((ftype)2.0/((ftype)(i+1)*(i+2))));
			}
			tgen.push(i+1,i+1,-std::sqrt((ftype)2.0*(ftype)(i+1)/(ftype)(i+2)));
			++igen;
		}


		// initialize temporary matrices and vectors to be used by memeber functions
		tmat=new ctype*[n];
		for(i=0; i<n; ++i) {
			tmat[i]=new ctype[n]();
		}
		tmat2=new ctype*[n];
		for(i=0; i<n; ++i) {
			tmat2[i]=new ctype[n]();
		}
		tvec=new ftype[ngen]();
	}

	~sun_algebra() {
		if(generator!=0) {
			delete[] generator;
		}
		ngen=0;
		for(int i=0; i<n; ++i) {
			delete[] tmat[i];
		}
		delete[] tmat;
		for(int i=0; i<n; ++i) {
			delete[] tmat2[i];
		}
		delete[] tmat2;
		delete[] tvec;
		n=0;
	}

	void get_alg_mat_ah(ftype* invec,ctype** outmat) {
		// set outmat[][] to be ii/2 * sum_{j=0}^{n*n-1} invec[i]*\lambda[i][][]
		// i.e. the anti-hermitian matrix representation of the Lie-algebra vector invec[]
		int i,j;
		for(i=0; i<n; ++i) {
			for(j=0; j<n; ++j) {
				outmat[i][j]=0;
			}
		}
		ctype nfc=ctype(0.0,0.5);
		ctype tiv;
		for(int igen=0; igen<ngen; ++igen) {
			tgen=generator+igen;
			telem0=tgen->elem;
			tiv=nfc*invec[igen];
			for(i=0; i<tgen->nelem; ++i) {
				telem=telem0+i;
				outmat[telem->ind1][telem->ind2]+=tiv*telem->val;
			}
		}
	}

	void get_alg_mat_h(ftype* invec,ctype** outmat) {
		// set outmat[][] to be 1/2 * sum_{j=0}^{n*n-1} invec[i]*\lambda[i][][]
		// i.e. the hermitian matrix representation of the Lie-algebra vector invec[]
		int i,j;
		for(i=0; i<n; ++i) {
			for(j=0; j<n; ++j) {
				outmat[i][j]=0;
			}
		}
		ftype nfc=0.5;
		ctype tiv;
		for(int igen=0; igen<ngen; ++igen) {
			tgen=generator+igen;
			telem0=tgen->elem;
			tiv=nfc*invec[igen];
			for(i=0; i<tgen->nelem; ++i) {
				telem=telem0+i;
				outmat[telem->ind1][telem->ind2]+=tiv*telem->val;
			}
		}
	}

	void get_grp_mat(ftype* invec,ctype** outmat) {
		// set outmat[][]=exp(ii/2 * sum_{j=0}^{n*n-1} invec[i]*\lambda[i][][])
		int i,j;
		for(i=0; i<n; ++i) {
			for(j=0; j<n; ++j) {
				tmat[i][j]=0;
			}
		}
		ctype nfc=ctype(0.0,0.5);
		ctype tiv;
		for(int igen=0; igen<ngen; ++igen) {
			tgen=generator+igen;
			telem0=tgen->elem;
			tiv=nfc*invec[igen];
			for(i=0; i<tgen->nelem; ++i) {
				telem=telem0+i;
				tmat[telem->ind1][telem->ind2]+=tiv*telem->val;
			}
		}
		ch(tmat,outmat);
	}

	void proj_ah(ctype** inmat,ftype* outvec) {
		// set outvec[] to be the vector of real parts of the coefficients obtained by projecting inmat[][]
		// on the anti-hermition Lie-algebra generators
		int i;
		ftype tiv;
		for(int igen=0; igen<ngen; ++igen) {
			tgen=generator+igen;
			telem0=tgen->elem;
			tiv=0;
			for(i=0; i<tgen->nelem; ++i) {
				telem=telem0+i;
				tiv+=std::imag(telem->val*inmat[telem->ind2][telem->ind1]);
			}
			outvec[igen]=tiv;
		}
	}

	void proj_ah(ctype** inmat,ctype** outmat) {
		// set outmat[] to be the anti-hermition Lie-algebra projection of inmat[][]
		int i,j;
		for(i=0; i<n; ++i) {
			for(j=0; j<n; ++j) {
				outmat[i][j]=0;
			}
		}
		ctype nfc=ctype(0.0,0.5);
		ftype tiv;
		ctype ttiv;
		for(int igen=0; igen<ngen; ++igen) {
			tgen=generator+igen;
			telem0=tgen->elem;
			tiv=0;
			for(i=0; i<tgen->nelem; ++i) {
				telem=telem0+i;
				tiv+=std::imag(telem->val*inmat[telem->ind2][telem->ind1]);
			}
			ttiv=tiv*nfc;
			for(i=0; i<tgen->nelem; ++i) {
				telem=telem0+i;
				outmat[telem->ind1][telem->ind2]+=ttiv*telem->val;
			}
		}
	}

	void proj_h(ctype** inmat,ftype* outvec) {
		// set outvec[] to be the vector of real parts of the coefficients obtained by projecting inmat[][]
		// on the hermition Lie-algebra generators
		int i;
		ftype tiv;
		for(int igen=0; igen<ngen; ++igen) {
			tgen=generator+igen;
			telem0=tgen->elem;
			tiv=0;
			for(i=0; i<tgen->nelem; ++i) {
				telem=telem0+i;
				tiv+=std::real(telem->val*inmat[telem->ind2][telem->ind1]);
			}
			outvec[igen]=tiv;
		}
	}

	void proj_h(ctype** inmat,ctype** outmat) {
		// set outmat[] to be the hermition Lie-algebra projection of inmat[][]

		int i,j;
		for(i=0; i<n; ++i) {
			for(j=0; j<n; ++j) {
				outmat[i][j]=0;
			}
		}
		ftype nfc=0.5;
		ftype tiv;
		for(int igen=0; igen<ngen; ++igen) {
			tgen=generator+igen;
			telem0=tgen->elem;
			tiv=0;
			for(i=0; i<tgen->nelem; ++i) {
				telem=telem0+i;
				tiv+=std::real(telem->val*inmat[telem->ind2][telem->ind1]);
			}
			tiv*=nfc;
			for(i=0; i<tgen->nelem; ++i) {
				telem=telem0+i;
				outmat[telem->ind1][telem->ind2]+=tiv*telem->val;
			}
		}
	}

	void get_log(ctype** inmat,ftype* outvec) {
		// computes the lie-algebra vector representing the matrix log of the SU(n) matrix inmat[][] (cf.arXiv:2404.07704)
		// (this implementation does the lie-algebra projection step via projection on the anti-hermitian lie-algebra generators) 
		int i;
		for(i=0; i<ngen; ++i) {
			outvec[i]=0; //make sure, outvec[] is zero-vector
		}
		ctype nfc=ctype(0.0,0.5);
		ftype tiv;
		ctype ttiv;
		int maxit=5*n; //limit the maximum number of iterations
		ftype tfprec=(ftype)10.0*_fprec*(ftype)n*n; //iteration will stop as soon as the 1-norm of the change done to outvec[]
											 //during the last iteration is smaller than tfprec times the 1-norm of outvec[] itself
		
		ch.matrix_copy(inmat,tmat); //start by setting tmat[][]=inmat[][]
		int it;
		ftype outvn=0;
		ftype tivn;
		for(it=0; it<maxit; ++it) {
			outvn=0;
			tivn=0;
			// add Lie-algebra projection of tmat[][] to outvec[]:
			for(int igen=0; igen<ngen; ++igen) {
				tgen=generator+igen;
				telem0=tgen->elem;
				tiv=0;
				for(i=0; i<tgen->nelem; ++i) {
					telem=telem0+i;
					tiv+=std::imag(telem->val*tmat[telem->ind2][telem->ind1]);
				}
				outvec[igen]+=tiv;
				outvn+=std::abs(outvec[igen]);
				tivn+=std::abs(tiv);
			}
			if(tivn<tfprec*outvn) {
				break;
			}

			// set tmat[][] to be the negative of the anti-hermitian matrix corresponding to outvec[]:
			ch.set_to_zero(tmat);
			for(int igen=0; igen<ngen; ++igen) {
				tgen=generator+igen;
				telem0=tgen->elem;
				ttiv=-nfc*outvec[igen];
				for(i=0; i<tgen->nelem; ++i) {
					telem=telem0+i;
					tmat[telem->ind1][telem->ind2]+=ttiv*telem->val;
				}
			}
			
			ch(tmat,tmat2); // set tmat2[][]=exp(tmat[][]); tmat2[][] is now an approximation to inverse(inmat)

			ch.matrix_mult_nn(inmat,tmat2,tmat); // set tmat[][]=inmat[][].tmat2[][]
		}
		//print out the number of iterations that were required and some additional info:
		//std::cout<<"iterations: "<<it<<" , outvn: "<<outvn<<"("<<outvn*tfprec<<") , tivn: "<<tivn<<" , _fprec: "<<_fprec<<std::endl;
	}


	int n;
	int ngen;
	sparse_mat<ctype>* generator;
	chexp<ctype> ch;
	ctype** tmat;
	ctype** tmat2;
	ftype* tvec;
	sparse_mat<ctype>* tgen;
	sa_entry<ctype>* telem0;
	sa_entry<ctype>* telem;
};
