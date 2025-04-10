/*
	header-only C++ implementation of iterative Cayley-Hamilton method 
	for computing matrix power series and their differentials (cf. arXiv:2404.07704)
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
#include<iostream>
#include<algorithm>
#include<complex>
#include<limits>
#include<cmath>
#include <iomanip>
#define FIXED_FLOAT(x) std::fixed<<std::setprecision(x)

typedef double ftype; // float type to be used
typedef long double lftype; // long version of ftype

//typedef float ftype; // float type to be used
//typedef double lftype; // long version of ftype

typedef std::complex<ftype> ctype; // complex type corresponding to ftype
typedef std::complex<lftype> lctype; // complex type corresponding to lftype

namespace type_helper {
	// if T is not a std::complex type, then val_type<T>::type will be T itself
	template<typename T>
	struct val_type {
		using type=T;
	};

	// if T is of type std::complex<fT>, then val_type<T>::type will be fT
	template<typename T>
	struct val_type<std::complex<T>> {
		using type=typename std::complex<T>::value_type;
	};

	// if no long version of type T is defined, then long_type<T>::type will be T itself
	template<typename T>
	struct long_type { typedef T type; };

	// if T is ftype, then long_type<T>::type will be lftype
	template<>
	struct long_type<ftype> { typedef lftype type; };

	// if T is ctype, then long_type<T>:: type will be lctype
	template<>
	struct long_type<ctype> { typedef lctype type; };
}

static const ftype _fprec=std::numeric_limits<ftype>::epsilon();

template<class T>
class cayley_hamilton {
	// template class providing an implementation of the iterative Cayley-Hamilton method for computing
	// matrix power series (cf. arXiv:2404.07704).

public:
	using fT=typename type_helper::val_type<T>::type; // underlying floating point type of type T
	using lT=typename type_helper::long_type<T>::type; // longer bit representation of type T (if defined) or type T itself

	cayley_hamilton(): n(0),pl(0),trpl(0),crpl(0),pal(0),al(0),mmax(0),nhl_max(0),tmat1(0),tmat2(0),kal(0) {
		//default constructor; will require a call to set_n() to set the size of the square matrices on which the class will operate

		opf=[](fT pref,int i) { return pref/(fT)i; }; // returns the i-th coeff. of the power seris, computed from the (i-1)-th coeff. "pref" and "i"
													        // default rule generates the coeffs. for the powerseries of exp(), i.e. 1/i!
	}

	cayley_hamilton(int tn,fT(*topf)(fT,int)=0) {
		// constructor with matrix size n=tn as argument and optional argument for function pointer or lambda that defines power series coeffients. 
		if(topf!=0) {
			// non-zero function pointer
			opf=topf;
		} else {
			opf=[](fT pref,int i) { return pref/(fT)i; };
		}
		if(tn>0) {
			// valid matrix size
			n=tn;
			new_matrix(tmat1);
			new_matrix(tmat2);
			new_matrix_array(pl,n+1);
			trpl=new T[n+1]();
			crpl=new T[n+1]();
			pal=new T[n]();
			al=new T[n]();
			kal=new T[n]();
			mmax=100*n;
			nhl_max=3;
		} else {
			// invalid matrix size --> default initialization
			n=0;
			tmat1=0;
			tmat2=0;
			pl=0;
			trpl=0;
			crpl=0;
			pal=0;
			al=0;
			kal=0;
			mmax=100;
			nhl_max=3;
		}
	}

	~cayley_hamilton() {
		// destructor
		if(tmat1!=0) {
			delete_matrix(tmat1);
		}
		if(tmat2!=0) {
			delete_matrix(tmat2);
		}
		if(pl!=0) {
			delete_matrix_array(pl,n+1);
		}
		if(trpl!=0) {
			delete[] trpl;
			trpl=0;
		}
		if(crpl!=0) {
			delete[] crpl;
			crpl=0;
		}
		if(pal!=0) {
			delete[] pal;
			pal=0;
		}
		if(al!=0) {
			delete[] al;
			al=0;
		}
		if(kal!=0) {
			delete[] kal;
			kal=0;
		}
		n=0;
	}

	void set_n(int tn) {
		//if an instance has been created with the trivial/empty constructor, 
		//this function has to be run afterwards to initialize the instance properly
		//and to prepare it to operate on matrices of size nxn with n=tn
		if(tn>0) {
			//tn is valid matrix size
			if(tn!=n) {
				//if tn differs from n, "resize" all the temporary arrays and matrices 
				if(tmat1!=0) {
					delete_matrix(tmat1);
				}
				if(tmat2!=0) {
					delete_matrix(tmat2);
				}
				if(pl!=0) {
					delete_matrix_array(pl,n+1);
				}
				if(trpl!=0) {
					delete[] trpl;
					trpl=0;
				}
				if(crpl!=0) {
					delete[] crpl;
					crpl=0;
				}
				if(pal!=0) {
					delete[] pal;
					pal=0;
				}
				if(al!=0) {
					delete[] al;
					al=0;
				}
				if(kal!=0) {
					delete[] kal;
					kal=0;
				}

				n=tn;

				new_matrix(tmat1);

				new_matrix(tmat2);

				new_matrix_array(pl,n+1);

				trpl=new T[n+1];

				crpl=new T[n+1];

				pal=new T[n];

				al=new T[n];

				kal=new T[n];
				
				mmax=100*n;
			}
		} else {
			// tn is not a valid matrix size --> free all memory
			if(tmat1!=0) {
				delete_matrix(tmat1);
			}
			if(tmat2!=0) {
				delete_matrix(tmat2);
			}
			if(pl!=0) {
				delete_matrix_array(pl,n+1);
			}
			if(trpl!=0) {
				delete[] trpl;
				trpl=0;
			}
			if(crpl!=0) {
				delete[] crpl;
				crpl=0;
			}
			if(pal!=0) {
				delete[] pal;
				pal=0;
			}
			if(al!=0) {
				delete[] al;
				al=0;
			}
			if(kal!=0) {
				delete[] kal;
				kal=0;
			}
			mmax=0;

			n=0;
		}
	}

	void set_opf(fT(*topf)(fT,int)) {
		//set opf to point to a user-defined function or lambda for generating the power series coefficients
		opf=topf;
	}

	int operator()(T** ain,T** aout,fT rescale=0) {
		// applies the power series defined by opf() to the matrix ain[][], using the Cayley-Hamilton recursion, 
		// and writes the result to aout[][]
		// if rescale is set to 1, then the computation will be performed with rescaled input matrix, which is
		// useful for matrices of large norm; if rescale isset to 0, no matrix rescaling will be performed. 
		int niter=0;
		if(n>0) {
			int i,j,k;
			fT sfac=(fT)1.0; //scaling factor
			// compute the 0-th to n-th matrix powers of ta[][] :
			//   the i-th matrix power of ta[][] is stored in pl[i][][]
			//   the trace of the i-th matrix power of ta[][] is stored in trpl[i]

			trpl[0]=n;
			if(rescale!=0) {
				// if matrix rescaling is used, set sfac to be the magnitude of the largest element of ain[][]
				// and initiate the computation of the matrix powers from ain[][]/sfac instead of ain[][]:
				// (note that since we compute also the matrix powers pl[] from the rescaled matrix,
				// the Cayley-Hamilton coefficient a_{k,j}, with k=0,1,...,\infty; j=0,...,n-1 
				// will need to be rescaled only by a factor  sfac^{k}, instead of sfac^{k-j})
				sfac=(fT)1.0/rescale;
				matrix_copy_scaled(ain,rescale,pl[1]);
			} else {
				matrix_copy(ain,pl[1]);
			}
			get_trace(pl[1],trpl[1]);
			for(i=2; i<n; ++i) {
				j=i/2;
				k=i%2;
				matrix_mult_nn(pl[j],pl[j+k],pl[i],trpl[i]);
			}
			j=n/2;
			k=n%2;
			matrix_mult_trace_nn(pl[j],pl[j+k],trpl[n]);

			// compute the characteristic polynomial crpl[] from the traced powers trpl[] :
			crpl[n]=1;
			for(j=1; j<=n; ++j) {
				crpl[n-j]=-trpl[j];
				for(i=1; i<j; ++i) {
					crpl[n-j]-=crpl[n-(j-i)]*trpl[i];
				}
				crpl[n-j]/=(fT)j;
			}

			// compute iteratively the n coefficients al[] so that the Cayley-Hamilton result
			// for the matrix power series is given by aout[][] = sum_i{al[i]*pl[i][][]} :

			// set initial values for the n entries in al[] and pal[] :
			fT wpf=1.0,twpf=1.0,wpfr=1.0,ttwpf; //leading coefficient of power series and of running sum used for convergence check
			for(i=0; i<n; ++i) {
				pal[i]=0;
				al[i]=wpf;
				wpf=opf(wpf,i+1); //compute (i+1)-th power series coefficent from the i-th coefficient, using the rule defined by opf()
				if(rescale!=0) {
					//if matrix rescaling is used, next power series term will need additional factor of sfac compared to current term.
					wpf*=sfac;
				}
				twpf+=wpf;
			}
			pal[n-1]=1.0;

			// next we iteratively add higher order power series terms to al[] till al[] stops changing
			// more precisely: the iteration will terminate after nhl_max consecutive iterations leave twpf unchanged
			T cho; // temporary variables for iterating
			int nhl=0; // counts the number of consecutive non-changing iterations  
			fT s,rs=1.0; // used for normalizing the pal[]
			for(j=n; j<mmax; ++j) {
				s=0;
				cho=pal[n-1]*rs;
				for(i=n-1; i>0; --i) {
					pal[i]=pal[i-1]*rs-cho*crpl[i];
					s+=std::norm(pal[i]);
					al[i]+=wpf*pal[i];
				}
				pal[0]=-cho*crpl[0];
				s+=std::norm(pal[0]);
				al[0]+=wpf*pal[0];

				wpf=opf(wpf,j+1); //compute (i+1)-th power series coefficent from i-th coefficient, using the rule defined by opf()

				if(rescale!=0) {
					//if matrix rescaling is used, next power series term will need additional factor of sfac compared to current term.
					wpf*=sfac;
				}

				s=std::sqrt(s);
				if(s>1.0) {
					// if s is bigger than 1, normalize pal[] by a factor rs=1.0/s in next itaration, and multiply wpf by s to compensate
					wpf*=s;
					wpfr=1.0;
					rs=1.0/s;
				} else {
					wpfr=s;
					rs=1.0;
				}

				ttwpf=twpf;
				twpf+=wpf*wpfr;
				if(ttwpf==twpf) {
					++nhl;
					if(nhl>=nhl_max) {
						//terminate iteration
						break;
					}
				} else {
					nhl=0;
				}

			}
			niter=j;

			// form output matrix by summing the 0-th to (n-1)-th matrix powers pl[] with corresponding weights al[] 
			set_to_identity_scaled(al[0],aout);
			for(j=1; j<n; ++j) {
				matrix_add_scaled(pl[j],al[j],aout);
			}
		}
		return niter;
	}

	int operator()(T** ain,T** aout,T** daout,fT rescale=0) {
		// applies the power series defined by opf() to the matrix ain[][], using the Cayley-Hamilton recursion, 
		// and writes the result to aout[][]; computes also the derivative of aout[][] in the direction of
		// daout[][] and overwrite daout[][] with the result :
		int niter=0;
		if(n>0) {
			// use unused entries of pl[] as storage for computation of differentaial components:
			T** kmats=pl[0];
			T** kh=pl[n];

			int i,j,k;
			fT sfac=(fT)1.0; //scaling factor
			// compute the 0-th to n-th matrix powers of ta[][] :
			//   the i-th matrix power of ta[][] is stored in pl[i][][]
			//   the trace of the i-th matrix power of ta[][] is stored in trpl[i]
			trpl[0]=n;
			if(rescale!=0) {
				// if matrix rescaling is used, set sfac to be the magnitude of the largest element of ain[][]
				// and initiate the computation of the matrix powers from ain[][]/sfac instead of ain[][]:
				// (note that since we compute also the matrix powers pl[] from the rescaled matrix,
				// the Cayley-Hamilton coefficient a_{k,j}, with k=0,1,...,\infty; j=0,...,n-1 
				// will need to be rescaled only by a factor  sfac^{k}, instead of sfac^{k-j})
				sfac=(fT)1.0/rescale;
				matrix_copy_scaled(ain,rescale,pl[1]);
			} else {
				matrix_copy(ain,pl[1]);
			}
			get_trace(pl[1],trpl[1]);
			for(i=2; i<n; ++i) {
				j=i/2;
				k=i%2;
				matrix_mult_nn(pl[j],pl[j+k],pl[i],trpl[i]);
			}
			j=n/2;
			k=n%2;
			matrix_mult_trace_nn(pl[j],pl[j+k],trpl[n]);

			// compute the characteristic polynomial crpl[] from the traced powers trpl[] :
			crpl[n]=1;
			for(j=1; j<=n; ++j) {
				crpl[n-j]=-trpl[j];
				for(i=1; i<j; ++i) {
					crpl[n-j]-=crpl[n-(j-i)]*trpl[i];
				}
				crpl[n-j]/=j;
			}

			// compute iteratively the n coefficients al[] so that the Cayley-Hamilton result
			// for the matrix power series is given by aout[][] = sum_i{al[i]*pl[i][][]} :

			// set initial values for the n entries in al[] and pal[] and the nxn entries in kmats[][] :
			fT wpf=1.0,twpf=1.0,wpfr=1.0,ttwpf; //leading coefficient of power series and of running sum used for convergence check
			set_to_zero(kmats);
			for(i=0; i<n; ++i) {
				pal[i]=0;
				al[i]=wpf;
				wpf=opf(wpf,i+1); //compute (i+1)-th power series coefficent from the i-th coefficient, using the rule defined by opf()

				k=i/2;
				for(j=i-k; j<=i; ++j) {
					kmats[i-j][j]=wpf;
				}

				if(rescale!=0) {
					//if matrix rescaling is used, next power series term will need additional factor of sfac compared to current term.
					wpf*=sfac;
				}
				twpf+=wpf;
			}
			pal[n-1]=1.0;
			// set initial values for the entries in kh[][] :
			set_to_zero(kh);
			k=(n-1)/2;
			for(i=n-1-k; i<n; ++i) {
				kh[n-1-i][i]=1.0;
			}


			// next we iteratively add higher order power series terms to al[] and kmats[][] till al[] stops changing
			// more precisely: the iteration will terminate after nhl_max consecutive iterations have not changed al[]	
			T cho; // temporary variables for iterating
			int nhl=0; // counts the number of consecutive non-changing iterations  
			fT s,rs=1.0; // used for normalizing the pal[]
			for(j=n; j<mmax; ++j) {
				s=0;
				cho=pal[n-1]*rs;
				for(i=n-1; i>0; --i) {
					pal[i]=pal[i-1]*rs-cho*crpl[i];
					s+=std::norm(pal[i]);
					al[i]+=wpf*pal[i];
				}
				pal[0]=-cho*crpl[0];
				s+=std::norm(pal[0]);
				al[0]+=wpf*pal[0];

				wpf=opf(wpf,j+1); //compute (i+1)-th power series coefficent from i-th coefficient, using the rule defined by opf()

				s=std::sqrt(s);
				if(s>1.0) {
					// if s is bigger than 1, normalize pal[] by a factor rs=1.0/s in next itaration, and multiply wpf by s to compensate
					wpf*=s;
					wpfr=1.0;
					rs=1.0/s;
				} else {
					wpfr=s;
					rs=1.0;
				}

				// add terms to kmats[][]:
				for(i=n-1; i>=0; --i) {
					cho=kh[i][n-1]*rs;
					for(k=n-1; k>i; --k) {
						kh[i][k]=kh[i][k-1]*rs-cho*crpl[k];
						kmats[i][k]+=wpf*kh[i][k];
					}
					if(i>0) {
						kh[i][i]=kh[i-1][i]*rs-cho*crpl[i];
					} else {
						kh[0][0]=pal[0]*rs-cho*crpl[0];
					}
					kmats[i][i]+=wpf*kh[i][i];
				}

				if(rescale!=0) {
					//if matrix rescaling is used, next power series term will need additional factor of sfac compared to current term.
					wpf*=sfac;
				}

				ttwpf=twpf;
				twpf+=wpf*wpfr;
				if(ttwpf==twpf) {
					++nhl;
					if(nhl>=nhl_max) {
						//terminate iteration
						break;
					}
				} else {
					nhl=0;
				}

			}
			niter=j;


			// form output matrix by summing the 0-th to (n-1)-th matrix powers pl[] with corresponding weights al[] 
			set_to_identity_scaled(al[0],aout);
			for(j=1; j<n; ++j) {
				matrix_add_scaled(pl[j],al[j],aout);
			}


			matrix_copy(daout,tmat1);

			//i=0, j=0:
			set_to_identity_scaled(kmats[0][0],kh);
			//i=0, j>0:
			for(j=1; j<n; ++j) {
				matrix_add_scaled(pl[j],kmats[0][j],kh);
			}
			matrix_mult_nn(tmat1,kh,daout);

			//i>0:
			for(i=1; i<n; ++i) {
				//j=0:
				set_to_identity_scaled(kmats[0][i],kh);
				//j>0:
				for(j=1; j<i; ++j) {
					matrix_add_scaled(pl[j],kmats[j][i],kh);
				}
				for(j=i; j<n; ++j) {
					matrix_add_scaled(pl[j],kmats[i][j],kh);
				}
				matrix_mult_nn(tmat1,kh,tmat2);
				matrix_mult_nn_add(pl[i],tmat2,daout);
			}

		}
		return niter;
	}

	int operator()(T** ain,T** aout,T**** daout,fT rescale=0) {
		// applies the power series defined by opf() to the matrix ain[][], using the Cayley-Hamilton recursion, 
		// and writes the result to aout[][]; computes also the derivative of aout[][] with respect to each of 
		// the nxn components of ain[][] and write the result to daout[][][][] (the first two indices define the component
		// with respect to which the derivative is taken and the last two indices enumerate the components of the matrix-
		// valued derivative).
		int niter=0;
		if(n>0) {
			// use unused entries of pl[] as storage for computation of differentaial components:
			T** kmats=pl[0];
			T** kh=pl[n];

			int i,j,k,l;
			fT sfac=1.0; //scaling factor
			// compute the 0-th to n-th matrix powers of ta[][] :
			//   the i-th matrix power of ta[][] is stored in pl[i][][]
			//   the trace of the i-th matrix power of ta[][] is stored in trpl[i]
			trpl[0]=n;
			if(rescale!=0) {
				// if matrix rescaling is used, set sfac to be the magnitude of the largest element of ain[][]
				// and initiate the computation of the matrix powers from ain[][]/sfac instead of ain[][]:
				// (note that since we compute also the matrix powers pl[] from the rescaled matrix,
				// the Cayley-Hamilton coefficient a_{k,j}, with k=0,1,...,\infty; j=0,...,n-1 
				// will need to be rescaled only by a factor  sfac^{k}, instead of sfac^{k-j})
				sfac=1.0/rescale;
				matrix_copy_scaled(ain,rescale,pl[1]);
			} else {
				matrix_copy(ain,pl[1]);
			}
			get_trace(pl[1],trpl[1]);
			for(i=2; i<n; ++i) {
				j=i/2;
				k=i%2;
				matrix_mult_nn(pl[j],pl[j+k],pl[i],trpl[i]);
			}
			j=n/2;
			k=n%2;
			matrix_mult_trace_nn(pl[j],pl[j+k],trpl[n]);

			// compute the characteristic polynomial crpl[] from the traced powers trpl[] :
			crpl[n]=1;
			for(j=1; j<=n; ++j) {
				crpl[n-j]=-trpl[j];
				for(i=1; i<j; ++i) {
					crpl[n-j]-=crpl[n-(j-i)]*trpl[i];
				}
				crpl[n-j]/=j;
			}

			// compute iteratively the n coefficients al[] so that the Cayley-Hamilton result
			// for the matrix power series is given by aout[][] = sum_i{al[i]*pl[i][][]} :

			// set initial values for the n entries in al[] and pal[] and the nxn entries in kmats[][] :
			fT wpf=1.0,twpf=1.0,wpfr=1.0,ttwpf; //leading coefficient of power series and of running sum used for convergence check
			set_to_zero(kmats);
			for(i=0; i<n; ++i) {
				pal[i]=0;
				al[i]=wpf;
				wpf=opf(wpf,i+1); //compute (i+1)-th power series coefficent from the i-th coefficient, using the rule defined by opf()

				k=i/2;
				for(j=i-k; j<=i; ++j) {
					kmats[i-j][j]=wpf;
				}

				if(rescale!=0) {
					//if matrix rescaling is used, next power series term will need additional factor of sfac compared to current term.
					wpf*=sfac;
				}
				twpf+=wpf;
			}
			pal[n-1]=1.0;
			// set initial values for the entries in kh[][] :
			set_to_zero(kh);
			k=(n-1)/2;
			for(i=n-1-k; i<n; ++i) {
				kh[n-1-i][i]=1.0;
			}


			// next we iteratively add higher order power series terms to al[] and kmats[][] till al[] stops changing
			// more precisely: the iteration will terminate after nhl_max consecutive iterations have not changed al[]	
			T cho; // temporary variables for iterating
			int nhl=0; // counts the number of consecutive non-changing iterations  
			fT s,rs=1.0; // used for normalizing the pal[]
			for(j=n; j<mmax; ++j) {
				s=0;
				cho=pal[n-1]*rs;
				for(i=n-1; i>0; --i) {
					pal[i]=pal[i-1]*rs-cho*crpl[i];
					s+=std::norm(pal[i]);
					al[i]+=wpf*pal[i];
				}
				pal[0]=-cho*crpl[0];
				s+=std::norm(pal[0]);
				al[0]+=wpf*pal[0];

				wpf=opf(wpf,j+1); //compute (i+1)-th power series coefficent from i-th coefficient, using the rule defined by opf()

				s=std::sqrt(s);
				if(s>1.0) {
					// if s is bigger than 1, normalize pal[] by a factor rs=1.0/s in next itaration, and multiply wpf by s to compensate
					wpf*=s;
					wpfr=1.0;
					rs=1.0/s;
				} else {
					wpfr=s;
					rs=1.0;
				}

				// add terms to kmats[][]:
				for(i=n-1; i>=0; --i) {
					cho=kh[i][n-1]*rs;
					for(k=n-1; k>i; --k) {
						kh[i][k]=kh[i][k-1]*rs-cho*crpl[k];
						kmats[i][k]+=wpf*kh[i][k];
					}
					if(i>0) {
						kh[i][i]=kh[i-1][i]*rs-cho*crpl[i];
					} else {
						kh[0][0]=pal[0]*rs-cho*crpl[0];
					}
					kmats[i][i]+=wpf*kh[i][i];
				}

				if(rescale!=0) {
					//if matrix rescaling is used, next power series term will need additional factor of sfac compared to current term.
					wpf*=sfac;
				}

				ttwpf=twpf;
				twpf+=wpf*wpfr;
				if(ttwpf==twpf) {
					++nhl;
					if(nhl>=nhl_max) {
						//terminate iteration
						break;
					}
				} else {
					nhl=0;
				}

			}
			niter=j;

			// form output matrix by summing the 0-th to (n-1)-th matrix powers pl[] with corresponding weights al[] 
			set_to_identity_scaled(al[0],aout);
			for(j=1; j<n; ++j) {
				matrix_add_scaled(pl[j],al[j],aout);
			}

			// form output for matrix of differentials daout[ic1][ic2][k][l] = daout[k][l]/dain[ic2][ic1] :
			// (note: in principle one could use symmetry domat[ic1][ic2].e(k, l) = domat[l][k].e(ic2, ic1),
			//  which would amount to  (i, j) <--> (j, i) with i = ic1 + n * ic2;  j = l + n * k; )

			//i=0, j=0:
			set_to_identity_scaled(kmats[0][0],kh);
			//i=0, j>0:
			for(j=1; j<n; ++j) {
				matrix_add_scaled(pl[j],kmats[0][j],kh);
			}
			int ic1,ic2;
			T** tdaout;
			for(ic1=0; ic1<n; ++ic1) {
				for(ic2=0; ic2<n; ++ic2) {
					tdaout=daout[ic1][ic2];
					set_to_zero(tdaout);
					for(k=0; k<n; ++k) {
						tdaout[k][ic1]+=kh[k][ic2];
					}
				}
			}

			//i>0:
			for(i=1; i<n; ++i) {
				//j=0:
				set_to_identity_scaled(kmats[0][i],kh);
				//j>0:
				for(j=1; j<i; ++j) {
					matrix_add_scaled(pl[j],kmats[j][i],kh);
				}
				for(j=i; j<n; ++j) {
					matrix_add_scaled(pl[j],kmats[i][j],kh);
				}
				for(ic1=0; ic1<n; ++ic1) {
					for(ic2=0; ic2<n; ++ic2) {
						tdaout=daout[ic1][ic2];
						for(k=0; k<n; ++k) {
							for(l=0; l<n; ++l) {
								tdaout[k][l]+=pl[i][ic1][l]*kh[k][ic2];
							}
						}

					}
				}

			}

		}
		return niter;
	}

	void get_r_k(T** ain,T* rout,T** kout,fT rescale=0) {
		// applies the power series defined by opf() to the matrix ain[][], using the Cayley-Hamilton recursion, 
		// and writes the Cayley-Hamilton coefficients to rout[]; computes also the decomposition matrix of the derivative
		// of rout[] with respect to the components of ain[][] and writes the result to kout[][] :
		if(n>0) {
			T* rl=rout;
			// use unused entries of pl[] as storage for computation of differentaial components:
			T** kmats=kout;
			T** kh=pl[n];

			int i,j,k;
			fT sfac=1.0; //scaling factor
			// compute the 0-th to n-th matrix powers of ta[][] :
			//   the i-th matrix power of ta[][] is stored in pl[i][][]
			//   the trace of the i-th matrix power of ta[][] is stored in trpl[i]
			trpl[0]=n;
			if(rescale!=0) {
				// if matrix rescaling is used, set sfac to be the magnitude of the largest element of ain[][]
				// and initiate the computation of the matrix powers from ain[][]/sfac instead of ain[][]:
				// (note that since we compute also the matrix powers pl[] from the rescaled matrix,
				// the Cayley-Hamilton coefficient a_{k,j}, with k=0,1,...,\infty; j=0,...,n-1 
				// will need to be rescaled only by a factor  sfac^{k}, instead of sfac^{k-j})
				sfac=1.0/rescale;
				matrix_copy_scaled(ain,rescale,pl[1]);
			} else {
				matrix_copy(ain,pl[1]);
			}
			get_trace(pl[1],trpl[1]);
			for(i=2; i<n; ++i) {
				j=i/2;
				k=i%2;
				matrix_mult_nn(pl[j],pl[j+k],pl[i],trpl[i]);
			}
			j=n/2;
			k=n%2;
			matrix_mult_trace_nn(pl[j],pl[j+k],trpl[n]);

			// compute the characteristic polynomial crpl[] from the traced powers trpl[] :
			crpl[n]=1;
			for(j=1; j<=n; ++j) {
				crpl[n-j]=-trpl[j];
				for(i=1; i<j; ++i) {
					crpl[n-j]-=crpl[n-(j-i)]*trpl[i];
				}
				crpl[n-j]/=j;
			}

			// compute iteratively the n coefficients al[] so that the Cayley-Hamilton result
			// for the matrix power series would be given by aout[][] = sum_i{al[i]*pl[i][][]} :

			// set initial values for the n entries in al[] and pal[] and the nxn entries in kmats[][] :
			fT wpf=1.0,twpf=1.0,wpfr=1.0,ttwpf; //leading coefficient of power series and of running sum used for convergence check
			set_to_zero(kmats);
			for(i=0; i<n; ++i) {
				pal[i]=0;
				rl[i]=wpf;
				wpf=opf(wpf,i+1); //compute (i+1)-th power series coefficent from the i-th coefficient, using the rule defined by opf()

				k=i/2;
				for(j=i-k; j<=i; ++j) {
					kmats[i-j][j]=wpf;
				}

				if(rescale!=0) {
					//if matrix rescaling is used, next power series term will need additional factor of sfac compared to current term.
					wpf*=sfac;
				}
				twpf+=wpf;
			}
			pal[n-1]=1.0;
			// set initial values for the entries in kh[][] :
			set_to_zero(kh);
			k=(n-1)/2;
			for(i=n-1-k; i<n; ++i) {
				kh[n-1-i][i]=1.0;
			}


			// next we iteratively add higher order power series terms to al[] and kmats[][] till al[] stops changing
			// more precisely: the iteration will terminate after nhl_max consecutive iterations have not changed al[]	
			T cho; // temporary variables for iterating
			int nhl=0; // counts the number of consecutive non-changing iterations  
			fT s,rs=1.0; // used for normalizing the pal[]
			for(j=n; j<mmax; ++j) {
				s=0;
				cho=pal[n-1]*rs;
				for(i=n-1; i>0; --i) {
					pal[i]=pal[i-1]*rs-cho*crpl[i];
					s+=std::norm(pal[i]);
					rl[i]+=wpf*pal[i];
				}
				pal[0]=-cho*crpl[0];
				s+=std::norm(pal[0]);
				rl[0]+=wpf*pal[0];

				wpf=opf(wpf,j+1); //compute (i+1)-th power series coefficent from i-th coefficient, using the rule defined by opf()

				s=std::sqrt(s);
				if(s>1.0) {
					// if s is bigger than 1, normalize pal[] by a factor rs=1.0/s in next itaration, and multiply wpf by s to compensate
					wpf*=s;
					wpfr=1.0;
					rs=1.0/s;
				} else {
					wpfr=s;
					rs=1.0;
				}

				// add terms to kmats[][]:
				for(i=n-1; i>=0; --i) {
					cho=kh[i][n-1]*rs;
					for(k=n-1; k>i; --k) {
						kh[i][k]=kh[i][k-1]*rs-cho*crpl[k];
						kmats[i][k]+=wpf*kh[i][k];
					}
					if(i>0) {
						kh[i][i]=kh[i-1][i]*rs-cho*crpl[i];
					} else {
						kh[0][0]=pal[0]*rs-cho*crpl[0];
					}
					kmats[i][i]+=wpf*kh[i][i];
				}

				if(rescale!=0) {
					//if matrix rescaling is used, next power series term will need additional factor of sfac compared to current term.
					wpf*=sfac;
				}

				ttwpf=twpf;
				twpf+=wpf*wpfr;
				if(ttwpf==twpf) {
					++nhl;
					if(nhl>=nhl_max) {
						//terminate iteration
						break;
					}
				} else {
					nhl=0;
				}

			}

		}
	}

	void get_r_dr(T** ain,T* rout,T*** drout,fT rescale=0) {
		// applies the power series defined by opf() to the matrix ain[][], using the Cayley-Hamilton recursion, 
		// and writes the Cayley-Hamilton coefficients to rout[]; computes also the derivatives of rout[]
		// with respect to the components of ain[][] and writes the resulting matrices to drout[][][] :
		if(n>0) {
			T* rl=rout;
			// use unused entries of pl[] as storage for computation of differentaial components:
			T** kmats=pl[0];
			T** kh=pl[n];

			int i,j,k;
			fT sfac=(fT)1.0; //scaling factor
			// compute the 0-th to n-th matrix powers of ta[][] :
			//   the i-th matrix power of ta[][] is stored in pl[i][][]
			//   the trace of the i-th matrix power of ta[][] is stored in trpl[i]
			trpl[0]=n;
			if(rescale!=0) {
				// if matrix rescaling is used, set sfac to be the magnitude of the largest element of ain[][]
				// and initiate the computation of the matrix powers from ain[][]/sfac instead of ain[][]:
				// (note that since we compute also the matrix powers pl[] from the rescaled matrix,
				// the Cayley-Hamilton coefficient a_{k,j}, with k=0,1,...,\infty; j=0,...,n-1 
				// will need to be rescaled only by a factor  sfac^{k}, instead of sfac^{k-j})
				sfac=(fT)1.0/rescale;
				matrix_copy_scaled(ain,rescale,pl[1]);
			} else {
				matrix_copy(ain,pl[1]);
			}
			get_trace(pl[1],trpl[1]);
			for(i=2; i<n; ++i) {
				j=i/2;
				k=i%2;
				matrix_mult_nn(pl[j],pl[j+k],pl[i],trpl[i]);
			}
			j=n/2;
			k=n%2;
			matrix_mult_trace_nn(pl[j],pl[j+k],trpl[n]);

			// compute the characteristic polynomial crpl[] from the traced powers trpl[] :
			crpl[n]=1;
			for(j=1; j<=n; ++j) {
				crpl[n-j]=-trpl[j];
				for(i=1; i<j; ++i) {
					crpl[n-j]-=crpl[n-(j-i)]*trpl[i];
				}
				crpl[n-j]/=j;
			}

			// compute iteratively the n coefficients al[] so that the Cayley-Hamilton result
			// for the matrix power series is given by aout[][] = sum_i{al[i]*pl[i][][]} :

			// set initial values for the n entries in al[] and pal[] and the nxn entries in kmats[][] :
			fT wpf=1.0,twpf=1.0,wpfr=1.0,ttwpf; //leading coefficient of power series and of running sum used for convergence check
			set_to_zero(kmats);
			for(i=0; i<n; ++i) {
				pal[i]=0;
				rl[i]=wpf;
				wpf=opf(wpf,i+1); //compute (i+1)-th power series coefficent from the i-th coefficient, using the rule defined by opf()

				k=i/2;
				for(j=i-k; j<=i; ++j) {
					kmats[i-j][j]=wpf;
				}

				if(rescale!=0) {
					//if matrix rescaling is used, next power series term will need additional factor of sfac compared to current term.
					wpf*=sfac;
				}
				twpf+=wpf;
			}
			pal[n-1]=1.0;
			// set initial values for the entries in kh[][] :
			set_to_zero(kh);
			k=(n-1)/2;
			for(i=n-1-k; i<n; ++i) {
				kh[n-1-i][i]=1.0;
			}


			// next we iteratively add higher order power series terms to al[] and kmats[][] till al[] stops changing
			// more precisely: the iteration will terminate after nhl_max consecutive iterations have not changed al[]	
			T cho; // temporary variables for iterating
			int nhl=0; // counts the number of consecutive non-changing iterations  
			fT s,rs=1.0; // used for normalizing the pal[]
			for(j=n; j<mmax; ++j) {
				s=0;
				cho=pal[n-1]*rs;
				for(i=n-1; i>0; --i) {
					pal[i]=pal[i-1]*rs-cho*crpl[i];
					s+=std::norm(pal[i]);
					rl[i]+=wpf*pal[i];
				}
				pal[0]=-cho*crpl[0];
				s+=std::norm(pal[0]);
				rl[0]+=wpf*pal[0];

				wpf=opf(wpf,j+1); //compute (i+1)-th power series coefficent from i-th coefficient, using the rule defined by opf()

				s=std::sqrt(s);
				if(s>1.0) {
					// if s is bigger than 1, normalize pal[] by a factor rs=1.0/s in next itaration, and multiply wpf by s to compensate
					wpf*=s;
					wpfr=1.0;
					rs=1.0/s;
				} else {
					wpfr=s;
					rs=1.0;
				}

				// add terms to kmats[][]:
				for(i=n-1; i>=0; --i) {
					cho=kh[i][n-1]*rs;
					for(k=n-1; k>i; --k) {
						kh[i][k]=kh[i][k-1]*rs-cho*crpl[k];
						kmats[i][k]+=wpf*kh[i][k];
					}
					if(i>0) {
						kh[i][i]=kh[i-1][i]*rs-cho*crpl[i];					
					} else {
						kh[0][0]=pal[0]*rs-cho*crpl[0];
					}
					kmats[i][i]+=wpf*kh[i][i];
				}

				if(rescale!=0) {
					//if matrix rescaling is used, next power series term will need additional factor of sfac compared to current term.
					wpf*=sfac;
				}

				ttwpf=twpf;
				twpf+=wpf*wpfr;
				if(ttwpf==twpf) {
					++nhl;
					if(nhl>=nhl_max) {
						//terminate iteration
						break;
					}
				} else {
					nhl=0;
				}

			}

			// compute drout[][][] from the kmats[][] matrix :

			// compute \tilde{R}-matrix:
			for(i=0; i<n; ++i) {
				tmat1[i][0]=kmats[i][n-1];
			}
			for(j=1; j<n; ++j) {
				cho=tmat1[n-1][j-1];
				for(i=n-1; i>0; --i) {
					tmat1[i][j]=tmat1[i-1][j-1]-cho*crpl[i];
				}
				tmat1[0][j]=-cho*crpl[0];
			}

			// compute R-matrix:
			for(j=0; j<n; ++j) {
				for(i=0; i<n; ++i) {
					tmat2[j][i]=tmat1[j][n-i-1];
					for(k=0; k<n-i-1; ++k) {
						tmat2[j][i]+=crpl[1+i+k]*tmat1[j][k];
					}
				}
			}
			// compute drout[]-matrices:
			for(k=0; k<n; ++k) {
				set_to_identity_scaled(tmat2[k][0],drout[k]);
				for(i=1; i<n; ++i) {
					matrix_mult_scalar_add(pl[i],tmat2[k][i],drout[k]);
				}
			}

		}
	}

	void ch_mult(T* tcrpl, T* ain1, T* ain2, T* aout) {
		// given the characteristic polynomial coefficients tcrpl[] for some matrix U, this routine
		// computes the Cayley-Hamilton decomposition coefficients of fout(U)=fin1(U).fin2(U) where
		// fin1(U) and fin2(U) are given by their Cayley-Hamilton decompositions ain1[] and ain2[]. 
		// the Cayley-Hamilton coefficients of fout(U) are written to aout[].
		int i,j;
		T cho=ain1[0];
		for(i=0; i<n; ++i) {
			kal[i]=ain1[i];
			pal[i]=ain2[i];
			aout[i]=cho*pal[i];
		}
		for(j=1; j<n; ++j) {
			cho=pal[n-1];
			for(i=n-1; i>0; --i) {
				pal[i]=pal[i-1]-cho*tcrpl[i];
				aout[i]+=kal[j]*pal[i];
			}
			pal[0]=-cho*tcrpl[0];
			aout[0]+=kal[j]*pal[0];
		}
	}

	// matrix operation utility functions:

	void set_to_zero(T** ta) {
		// set ta[][] to the zero matrix
		int ic1,ic2;
		T* l01;
		for(ic1=0; ic1<n; ++ic1) {
			l01=ta[ic1];
			for(ic2=0; ic2<n; ++ic2) {
				l01[ic2]=0;
			}
		}
	}

	void set_to_identity(T** ta) {
		// set ta[][] to the identiy matrix
		int ic1,ic2;
		T* l01;
		for(ic1=0; ic1<n; ++ic1) {
			l01=ta[ic1];
			for(ic2=0; ic2<n; ++ic2) {
				l01[ic2]=0;
			}
			l01[ic1]=1.0;
		}
	}

	template<typename sT>
	void set_to_identity_scaled(const sT& scalef, T** ta) {
		// set ta[][] to the identiy matrix
		int ic1,ic2;
		T* l01;
		for(ic1=0; ic1<n; ++ic1) {
			l01=ta[ic1];
			for(ic2=0; ic2<n; ++ic2) {
				l01[ic2]=0;
			}
			l01[ic1]=scalef;
		}
	}

	void matrix_copy(T** lin1,T** lout) {
		// set lout[][]=lin1[][]
		int ic1,ic2;
		T* lin10;
		T* lout0;
		for(ic1=0; ic1<n; ++ic1) {
			lout0=lout[ic1];
			lin10=lin1[ic1];
			for(ic2=0; ic2<n; ++ic2) {
				lout0[ic2]=lin10[ic2];
			}
		}
	}

	template<typename sT>
	void matrix_copy_scaled(T** lin1,const sT& scalef,T** lout) {
		// set lout[][]=scalef*lin1[][]
		int ic1,ic2;
		T* lin10;
		T* lout0;
		for(ic1=0; ic1<n; ++ic1) {
			lout0=lout[ic1];
			lin10=lin1[ic1];
			for(ic2=0; ic2<n; ++ic2) {
				lout0[ic2]=lin10[ic2]*scalef;
			}
		}
	}

	void matrix_copy_a(T** lin1,T** lout) {
		// set lout[][]=conjugate_transpose(lin1[][])
		int ic1,ic2;
		T* lin10;
		for(ic1=0; ic1<n; ++ic1) {
			lin10=lin1[ic1];
			for(ic2=0; ic2<n; ++ic2) {
				lout[ic2][ic1]=std::conj(lin10[ic2]);
			}
		}
	}

	template<typename sT>
	void matrix_copy_a_scaled(T** lin1,const sT& scalef,T** lout) {
		// set lout[][]=scalef*conjugate_transpose(lin1[][])
		int ic1,ic2;
		T* lin10;
		for(ic1=0; ic1<n; ++ic1) {
			lin10=lin1[ic1];
			for(ic2=0; ic2<n; ++ic2) {
				lout[ic2][ic1]=scalef*std::conj(lin10[ic2]);
			}
		}
	}

	template<typename sT>
	void matrix_add(const sT& scalef,T** lout) {
		for(int ic1=0; ic1<n; ++ic1) {
			lout[ic1][ic1]+=scalef;
		}
	}

	template<typename sT>
	void matrix_sub(const sT& scalef,T** lout) {
		for(int ic1=0; ic1<n; ++ic1) {
			lout[ic1][ic1]-=scalef;
		}
	}

	void matrix_add(T** lin1,T** lout) {
		// add lin1[][] to lout[][]
		int ic1,ic2;
		T* lin10;
		T* lout0;
		for(ic1=0; ic1<n; ++ic1) {
			lout0=lout[ic1];
			lin10=lin1[ic1];
			for(ic2=0; ic2<n; ++ic2) {
				lout0[ic2]+=lin10[ic2];
			}
		}
	}

	void matrix_add(T** lin1,T** lin2,T** lout) {
		// set lout[][] = lin1[][] + lin2[][] 
		int ic1,ic2;
		T* lin10;
		T* lin20;
		T* lout0;
		for(ic1=0; ic1<n; ++ic1) {
			lout0=lout[ic1];
			lin10=lin1[ic1];
			lin20=lin2[ic1];
			for(ic2=0; ic2<n; ++ic2) {
				lout0[ic2]=lin10[ic2]+lin20[ic2];
			}
		}
	}

	void matrix_add(T** lin1,T** lout, T& trout) {
		// add lin1[][] to lout[][] and set trout=trace(lout[][])
		int ic1,ic2;
		T* lin10;
		T* lout0;
		trout=0;
		for(ic1=0; ic1<n; ++ic1) {
			lout0=lout[ic1];
			lin10=lin1[ic1];
			for(ic2=0; ic2<n; ++ic2) {
				lout0[ic2]+=lin10[ic2];
			}
			trout+=lout0[ic1];
		}
	}

	void matrix_sub(T** lin1,T** lout) {
		// subtract lin1[][] from lout[][]
		int ic1,ic2;
		T* lin10;
		T* lout0;
		for(ic1=0; ic1<n; ++ic1) {
			lout0=lout[ic1];
			lin10=lin1[ic1];
			for(ic2=0; ic2<n; ++ic2) {
				lout0[ic2]-=lin10[ic2];
			}
		}
	}

	void matrix_sub(T** lin1,T** lin2,T** lout) {
		// set lout[][] = lin1[][] - lin2[][] 
		int ic1,ic2;
		T* lin10;
		T* lin20;
		T* lout0;
		for(ic1=0; ic1<n; ++ic1) {
			lout0=lout[ic1];
			lin10=lin1[ic1];
			lin20=lin2[ic1];
			for(ic2=0; ic2<n; ++ic2) {
				lout0[ic2]=lin10[ic2]-lin20[ic2];
			}
		}
	}

	void matrix_sub(T** lin1,T** lout,T& trout) {
		// subtract lin1[][] from lout[][] and set trout=trace(lout[][])
		int ic1,ic2;
		T* lin10;
		T* lout0;
		trout=0;
		for(ic1=0; ic1<n; ++ic1) {
			lout0=lout[ic1];
			lin10=lin1[ic1];
			for(ic2=0; ic2<n; ++ic2) {
				lout0[ic2]-=lin10[ic2];
			}
			trout+=lout0[ic1];
		}
	}

	template<typename sT>
	void matrix_add_scaled(T** lin1,const sT& scalef,T** lout) {
		// add scalef*lin1[][] to lout[][]
		int ic1,ic2;
		T* lin10;
		T* lout0;
		for(ic1=0; ic1<n; ++ic1) {
			lout0=lout[ic1];
			lin10=lin1[ic1];
			for(ic2=0; ic2<n; ++ic2) {
				lout0[ic2]+=scalef*lin10[ic2];
			}
		}
	}

	template<typename sT>
	void matrix_add_scaled(T** lin1,const sT& scalef,T** lout, T& trout) {
		// add scalef*lin1[][] to lout[][] and set trout=trace(lout[][])
		int ic1,ic2;
		T* lin10;
		T* lout0;
		trout=0;
		for(ic1=0; ic1<n; ++ic1) {
			lout0=lout[ic1];
			lin10=lin1[ic1];
			for(ic2=0; ic2<n; ++ic2) {
				lout0[ic2]+=scalef*lin10[ic2];
			}
			trout+=lout0[ic1];
		}
	}

	void matrix_mult_trace_nn(T** lin1,T** lin2,T& trout) {
		// set trout = trace(lin1[][].lin2[][])
		int ic1,ic2;
		T* lin10;
		trout=0;
		for(ic1=0; ic1<n; ++ic1) {
			lin10=lin1[ic1];
			for(ic2=0; ic2<n; ++ic2) {
				trout+=lin10[ic2]*lin2[ic2][ic1];
			}
		}
	}

	void matrix_mult_trace_na(T** lin1,T** lin2,T& trout) {
		// set trout = trace(lin1[][].lin2[][].conjugate_transpose())
		int ic1,ic2;
		T* lin10;
		T* lin20;
		trout=0;
		for(ic1=0; ic1<n; ++ic1) {
			lin10=lin1[ic1];
			lin20=lin2[ic1];
			for(ic2=0; ic2<n; ++ic2) {
				trout+=lin10[ic2]*std::conj(lin20[ic2]);
			}
		}
	}

	void matrix_mult_trace_an(T** lin1,T** lin2,T& trout) {
		// set trout = trace(lin1[][].lin2[][].conjugate_transpose())
		int ic1,ic2;
		T* lin10;
		T* lin20;
		trout=0;
		for(ic1=0; ic1<n; ++ic1) {
			lin10=lin1[ic1];
			lin20=lin2[ic1];
			for(ic2=0; ic2<n; ++ic2) {
				trout+=std::conj(lin10[ic2])*lin20[ic2];
			}
		}
	}

	void matrix_mult_trace_aa(T** lin1,T** lin2,T& trout) {
		// set trout = trace(lin1[][].lin2[][]).conj()
		int ic1,ic2;
		T* lin10;
		trout=0;
		for(ic1=0; ic1<n; ++ic1) {
			lin10=lin1[ic1];
			for(ic2=0; ic2<n; ++ic2) {
				trout+=lin10[ic2]*lin2[ic2][ic1];
			}
		}
		trout=std::conj(trout);
	}

	void matrix_mult_nn(T** lin1,T** lin2,T** lout) {
		// set lout[][] = lin1[][].lin2[][] (matrix product)
		int ic1,ic2,ic3;
		T* lin10;
		T* lout0;
		T tout;
		for(ic1=0; ic1<n; ++ic1) {
			lout0=lout[ic1];
			lin10=lin1[ic1];
			for(ic2=0; ic2<n; ++ic2) {
				tout=0;
				for(ic3=0; ic3<n; ++ic3) {
					tout+=lin10[ic3]*lin2[ic3][ic2];
				}
				lout0[ic2]=tout;
			}
		}
	}

	void matrix_mult_nn(T** lin1,T** lin2,T** lout,T& trout) {
		// set lout[][] = lin1[][].lin2[][] and trout=trace(lout[][])
		int ic1,ic2,ic3;
		T* lin10;
		T* lout0;
		T tout;
		trout=0;
		for(ic1=0; ic1<n; ++ic1) {
			lout0=lout[ic1];
			lin10=lin1[ic1];
			for(ic2=0; ic2<n; ++ic2) {
				tout=0;
				for(ic3=0; ic3<n; ++ic3) {
					tout+=lin10[ic3]*lin2[ic3][ic2];
				}
				lout0[ic2]=tout;
			}
			trout+=lout0[ic1];
		}
	}

	void matrix_mult_na(T** lin1,T** lin2,T** lout) {
		// set lout[][] = lin1[][].conjugate_transpose(lin2[][])
		int ic1,ic2,ic3;
		T* lin10;
		T* lin20;
		T* lout0;
		T tout;
		for(ic1=0; ic1<n; ++ic1) {
			lout0=lout[ic1];
			lin10=lin1[ic1];
			for(ic2=0; ic2<n; ++ic2) {
				lin20=lin2[ic2];
				tout=0;
				for(ic3=0; ic3<n; ++ic3) {
					tout+=lin10[ic3]*std::conj(lin20[ic3]);
				}
				lout0[ic2]=tout;
			}
		}
	}

	void matrix_mult_na(T** lin1,T** lin2,T** lout,T& trout) {
		// set lout[][] = lin1[][].conjugate_transpose(lin2[][]) and trout=trace(lout[][])
		int ic1,ic2,ic3;
		T* lin10;
		T* lin20;
		T* lout0;
		T tout;
		trout=0;
		for(ic1=0; ic1<n; ++ic1) {
			lout0=lout[ic1];
			lin10=lin1[ic1];
			for(ic2=0; ic2<n; ++ic2) {
				lin20=lin2[ic2];
				tout=0;
				for(ic3=0; ic3<n; ++ic3) {
					tout+=lin10[ic3]*std::conj(lin20[ic3]);
				}
				lout0[ic2]=tout;
			}
			trout+=lout0[ic1];
		}
	}

	void matrix_mult_an(T** lin1,T** lin2,T** lout) {
		// set lout[][] = conjugate_transpose(lin1[][]).lin2[][]
		int ic1,ic2,ic3;
		T* lout0;
		T tout;
		for(ic1=0; ic1<n; ++ic1) {
			lout0=lout[ic1];
			for(ic2=0; ic2<n; ++ic2) {
				tout=0;
				for(ic3=0; ic3<n; ++ic3) {
					tout+=std::conj(lin1[ic3][ic1])*lin2[ic3][ic2];
				}
				lout0[ic2]=tout;
			}
		}
	}

	void matrix_mult_an(T** lin1,T** lin2,T** lout,T& trout) {
		// set lout[][] = conjugate_transpose(lin1[][]).lin2[][] and trout=trace(lout[][])
		int ic1,ic2,ic3;
		T* lout0;
		T tout;
		trout=0;
		for(ic1=0; ic1<n; ++ic1) {
			lout0=lout[ic1];
			for(ic2=0; ic2<n; ++ic2) {
				tout=0;
				for(ic3=0; ic3<n; ++ic3) {
					tout+=std::conj(lin1[ic3][ic1])*lin2[ic3][ic2];
				}
				lout0[ic2]=tout;
			}
			trout+=lout0[ic1];
		}
	}

	template<typename sT>
	void matrix_mult_scalar(T** lin1,const sT& scalef,T** lout) {
		// add lin1[][]*scalef to lout[][]
		int ic1,ic2;
		T* lin10;
		T* lout0;
		for(ic1=0; ic1<n; ++ic1) {
			lout0=lout[ic1];
			lin10=lin1[ic1];
			for(ic2=0; ic2<n; ++ic2) {
				lout0[ic2]=lin10[ic2]*(T)scalef;
			}
		}
	}

	template<typename sT>
	void matrix_mult_scalar_add(T** lin1,const sT& scalef,T** lout) {
		// add lin1[][]*scalef to lout[][]
		int ic1,ic2;
		T* lin10;
		T* lout0;
		for(ic1=0; ic1<n; ++ic1) {
			lout0=lout[ic1];
			lin10=lin1[ic1];
			for(ic2=0; ic2<n; ++ic2) {
				lout0[ic2]+=lin10[ic2]*scalef;
			}
		}
	}

	template<typename sT>
	void matrix_mult_scalar_sub(T** lin1,const sT& scalef,T** lout) {
		// add lin1[][]*scalef to lout[][]
		int ic1,ic2;
		T* lin10;
		T* lout0;
		for(ic1=0; ic1<n; ++ic1) {
			lout0=lout[ic1];
			lin10=lin1[ic1];
			for(ic2=0; ic2<n; ++ic2) {
				lout0[ic2]-=lin10[ic2]*scalef;
			}
		}
	}

	void matrix_mult_nn_add(T** lin1,T** lin2,T** lout) {
		// add lin1[][].lin2[][] to lout[][]
		int ic1,ic2,ic3;
		T* lin10;
		T* lout0;
		T tout;
		for(ic1=0; ic1<n; ++ic1) {
			lout0=lout[ic1];
			lin10=lin1[ic1];
			for(ic2=0; ic2<n; ++ic2) {
				tout=0;
				for(ic3=0; ic3<n; ++ic3) {
					tout+=lin10[ic3]*lin2[ic3][ic2];
				}
				lout0[ic2]+=tout;
			}
		}
	}

	void matrix_mult_nn_add(T** lin1,T** lin2,T** lout,T& trout) {
		// add lin1[][].lin2[][] to lout[][] and set trout=trace(lout[][])
		int ic1,ic2,ic3;
		T* lin10;
		T* lout0;
		T tout;
		trout=0;
		for(ic1=0; ic1<n; ++ic1) {
			lout0=lout[ic1];
			lin10=lin1[ic1];
			for(ic2=0; ic2<n; ++ic2) {
				tout=0;
				for(ic3=0; ic3<n; ++ic3) {
					tout+=lin10[ic3]*lin2[ic3][ic2];
				}
				lout0[ic2]+=tout;
			}
			trout+=lout0[ic1];
		}
	}

	void matrix_mult_na_add(T** lin1,T** lin2,T** lout) {
		// add lin1[][].conjugate_transpose(lin2[][]) to lout[][]
		int ic1,ic2,ic3;
		T* lin10;
		T* lin20;
		T* lout0;
		T tout;
		for(ic1=0; ic1<n; ++ic1) {
			lout0=lout[ic1];
			lin10=lin1[ic1];
			for(ic2=0; ic2<n; ++ic2) {
				lin20=lin2[ic2];
				tout=0;
				for(ic3=0; ic3<n; ++ic3) {
					tout+=lin10[ic3]*std::conj(lin20[ic3]);
				}
				lout0[ic2]+=tout;
			}
		}
	}

	void matrix_mult_na_add(T** lin1,T** lin2,T** lout,T& trout) {
		// add lin1[][].conjugate_transpose(lin2[][]) to lout[][] and set trout=trace(lout[][])
		int ic1,ic2,ic3;
		T* lin10;
		T* lin20;
		T* lout0;
		T tout;
		trout=0;
		for(ic1=0; ic1<n; ++ic1) {
			lout0=lout[ic1];
			lin10=lin1[ic1];
			for(ic2=0; ic2<n; ++ic2) {
				lin20=lin2[ic2];
				tout=0;
				for(ic3=0; ic3<n; ++ic3) {
					tout+=lin10[ic3]*std::conj(lin20[ic3]);
				}
				lout0[ic2]+=tout;
			}
			trout+=lout0[ic1];
		}
	}

	void matrix_mult_an_add(T** lin1,T** lin2,T** lout) {
		// add conjugate_transpose(lin1[][]).lin2[][] to lout[][]
		int ic1,ic2,ic3;
		T* lout0;
		T tout;
		for(ic1=0; ic1<n; ++ic1) {
			lout0=lout[ic1];
			for(ic2=0; ic2<n; ++ic2) {
				tout=0;
				for(ic3=0; ic3<n; ++ic3) {
					tout+=std::conj(lin1[ic3][ic1])*lin2[ic3][ic2];
				}
				lout0[ic2]+=tout;
			}
		}
	}

	void matrix_mult_an_add(T** lin1,T** lin2,T** lout,T& trout) {
		// add conjugate_transpose(lin1[][]).lin2[][] to lout[][] and set trout=trace(lout[][])
		int ic1,ic2,ic3;
		T* lout0;
		T tout;
		trout=0;
		for(ic1=0; ic1<n; ++ic1) {
			lout0=lout[ic1];
			for(ic2=0; ic2<n; ++ic2) {
				tout=0;
				for(ic3=0; ic3<n; ++ic3) {
					tout+=std::conj(lin1[ic3][ic1])*lin2[ic3][ic2];
				}
				lout0[ic2]+=tout;
			}
			trout+=lout0[ic1];
		}
	}

	void matrix_mult_nn_sub(T** lin1,T** lin2,T** lout) {
		// subtract lin1[][].lin2[][] from lout[][]
		int ic1,ic2,ic3;
		T* lin10;
		T* lout0;
		T tout;
		for(ic1=0; ic1<n; ++ic1) {
			lout0=lout[ic1];
			lin10=lin1[ic1];
			for(ic2=0; ic2<n; ++ic2) {
				tout=0;
				for(ic3=0; ic3<n; ++ic3) {
					tout+=lin10[ic3]*lin2[ic3][ic2];
				}
				lout0[ic2]-=tout;
			}
		}
	}

	void matrix_mult_nn_sub(T** lin1,T** lin2,T** lout,T& trout) {
		// subtract lin1[][].lin2[][] from lout[][] and set trout=trace(lout[][])
		int ic1,ic2,ic3;
		T* lin10;
		T* lout0;
		T tout;
		trout=0;
		for(ic1=0; ic1<n; ++ic1) {
			lout0=lout[ic1];
			lin10=lin1[ic1];
			for(ic2=0; ic2<n; ++ic2) {
				tout=0;
				for(ic3=0; ic3<n; ++ic3) {
					tout+=lin10[ic3]*lin2[ic3][ic2];
				}
				lout0[ic2]-=tout;
			}
			trout+=lout0[ic1];
		}
	}

	void matrix_mult_na_sub(T** lin1,T** lin2,T** lout) {
		// subtract lin1[][].conjugate_transpose(lin2[][]) from lout[][]
		int ic1,ic2,ic3;
		T* lin10;
		T* lin20;
		T* lout0;
		T tout;
		for(ic1=0; ic1<n; ++ic1) {
			lout0=lout[ic1];
			lin10=lin1[ic1];
			for(ic2=0; ic2<n; ++ic2) {
				lin20=lin2[ic2];
				tout=0;
				for(ic3=0; ic3<n; ++ic3) {
					tout+=lin10[ic3]*std::conj(lin20[ic3]);
				}
				lout0[ic2]-=tout;
			}
		}
	}

	void matrix_mult_na_sub(T** lin1,T** lin2,T** lout,T& trout) {
		// subtract lin1[][].conjugate_transpose(lin2[][]) from lout[][] and set trout=trace(lout[][])
		int ic1,ic2,ic3;
		T* lin10;
		T* lin20;
		T* lout0;
		T tout;
		trout=0;
		for(ic1=0; ic1<n; ++ic1) {
			lout0=lout[ic1];
			lin10=lin1[ic1];
			for(ic2=0; ic2<n; ++ic2) {
				lin20=lin2[ic2];
				tout=0;
				for(ic3=0; ic3<n; ++ic3) {
					tout+=lin10[ic3]*std::conj(lin20[ic3]);
				}
				lout0[ic2]-=tout;
			}
			trout+=lout0[ic1];
		}
	}

	void matrix_mult_an_sub(T** lin1,T** lin2,T** lout) {
		// subtract conjugate_transpose(lin1[][]).lin2[][] from lout[][]
		int ic1,ic2,ic3;
		T* lout0;
		T tout;
		for(ic1=0; ic1<n; ++ic1) {
			lout0=lout[ic1];
			for(ic2=0; ic2<n; ++ic2) {
				tout=0;
				for(ic3=0; ic3<n; ++ic3) {
					tout+=std::conj(lin1[ic3][ic1])*lin2[ic3][ic2];
				}
				lout0[ic2]-=tout;
			}
		}
	}

	void matrix_mult_an_sub(T** lin1,T** lin2,T** lout,T& trout) {
		// subtract conjugate_transpose(lin1[][]).lin2[][] from lout[][] and set trout=trace(lout[][])
		int ic1,ic2,ic3;
		T* lout0;
		T tout;
		trout=0;
		for(ic1=0; ic1<n; ++ic1) {
			lout0=lout[ic1];
			for(ic2=0; ic2<n; ++ic2) {
				tout=0;
				for(ic3=0; ic3<n; ++ic3) {
					tout+=std::conj(lin1[ic3][ic1])*lin2[ic3][ic2];
				}
				lout0[ic2]-=tout;
			}
			trout+=lout0[ic1];
		}
	}

	void get_trace(T** lin1,T& rout) {
		// set rout=trace(lin1[][])
		rout=0;
		for(int ic1=0; ic1<n; ++ic1) {
			rout+=lin1[ic1][ic1];
		}
	}

	void get_frobenius_norm(T** lin1,fT& rout) {
		// set rout=frobenius_norm(lin1[][])
		rout=0;
		int ic1,ic2;
		for(ic1=0; ic1<n; ++ic1) {
			for(ic2=0; ic2<n; ++ic2) {
				rout+=std::norm(lin1[ic1][ic2]);
			}
		}
		rout=std::sqrt(rout);
	}

	void get_one_norm(T** lin1,fT& rout) {
		// set rout=frobenius_norm(lin1[][])
		rout=0;
		int ic1,ic2;
		for(ic1=0; ic1<n; ++ic1) {
			for(ic2=0; ic2<n; ++ic2) {
				rout+=std::abs(lin1[ic1][ic2]);
			}
		}
		rout=std::sqrt(rout);
	}

	void get_l1_norm(T** lin1,fT& rout) {
		// set rout=l1_norm(lin1[][])
		rout=0;
		int ic1,ic2;
		fT cols;
		for(ic2=0; ic2<n; ++ic2) {
			cols=0;
			for(ic1=0; ic1<n; ++ic1) {
				cols+=std::abs(lin1[ic1][ic2]);
			}
			if(cols>rout) {
				rout=cols;
			}
		}
	}

	void get_li_norm(T** lin1,fT& rout) {
		// set rout=li_norm(lin1[][])
		rout=0;
		int ic1,ic2;
		fT rows;
		for(ic1=0; ic1<n; ++ic1) {
			rows=0;
			for(ic2=0; ic2<n; ++ic2) {
				rows+=std::abs(lin1[ic1][ic2]);
			}
			if(rows>rout) {
				rout=rows;
			}
		}
	}

	void new_matrix(T**& ta,int init=-1) {
		// allocate memory for a nxn matrix and set ta to point to it
		ta=new T*[n];
		if(init==-1) {
			// do not initialize
			for(int ic1=0; ic1<n; ++ic1) {
				ta[ic1]=new T[n];
			}
		} else if(init==0) {
			// initialize with zero-matrix
			for(int ic1=0; ic1<n; ++ic1) {
				ta[ic1]=new T[n]();
			}
		} else {
			// initialize with identity matrix
			for(int ic1=0; ic1<n; ++ic1) {
				ta[ic1]=new T[n]();
				ta[ic1][ic1]=1;
			}
		}
	}

	void delete_matrix(T**& ta) {
		// deallocate storage pointed at by ta and set ta=0
		for(int ic1=0; ic1<n; ++ic1) {
			delete[] ta[ic1];
		}
		delete[] ta;
		ta=0;
	}

	void new_matrix_array(T***& tnna,int len,int init=-1) {
		// allocate memory for an array of nxn matrices of length "len"
		tnna=new T**[len];
		if(init==-1) {
			// do not initialize the matrices
			for(int i=0; i<len; ++i) {
				tnna[i]=new T*[n];
				for(int ic1=0; ic1<n; ++ic1) {
					tnna[i][ic1]=new T[n];
				}
			}
		} else if(init==0) {
			// initialize matrices to be zero matrices
			int ic1;
			for(int i=0; i<len; ++i) {
				tnna[i]=new T*[n];
				for(ic1=0; ic1<n; ++ic1) {
					tnna[i][ic1]=new T[n]();
				}
			}
		} else {
			// initialize matrices to be identity matrices
			int ic1;
			for(int i=0; i<len; ++i) {
				tnna[i]=new T*[n];
				for(ic1=0; ic1<n; ++ic1) {
					tnna[i][ic1]=new T[n]();
					tnna[i][ic1][ic1]=1;
				}
			}
		}
	}

	void delete_matrix_array(T***& tnn,int len) {
		// deallocate array of length "len" of nxn matrices pointed at by tnn and set tnn to zero
		for(int i=0; i<len; ++i) {
			for(int ic1=0; ic1<n; ++ic1) {
				delete[] tnn[i][ic1];
			}
			delete[] tnn[i];
		}
		delete[] tnn;
		tnn=0;
	}

	void print_matrix(T** ta, int dprec=6) {
		int totw=dprec+6;
		std::cout<<FIXED_FLOAT(dprec);
		
		for(int i=0; i<n; ++i) {
			for(int j=0; j<n; ++j) {
				std::cout.width(totw); std::cout<<ta[i][j]<<" ";
			}
			std::cout<<std::endl;
		}
	}

	void print_matrix_e(T** ta,int dprec=6) {
		std::cout<<FIXED_FLOAT(dprec);
		for(int i=0; i<n; ++i) {
			for(int j=0; j<n; ++j) {
				printf("% 0.*e ",dprec,ta[i][j]);
			}
			std::cout<<std::endl;
		}
	}

protected:
	int n;
	int mmax;
	int nhl_max;
	T*** pl;
	T* trpl;
	T* crpl;
	T* pal;
	T* al;
	T** tmat1;
	T** tmat2;
	T* kal;
private:
	fT(*opf)(fT,int);
};

template<class T>
class chexp : public cayley_hamilton<T> {
// template class providing an implementation of matrix exponentiation using iterative Cayley-Hamilton method 
// with scaling and squaring. (cf. arXiv:2404.07704).
public:
	using fT=typename cayley_hamilton<T>::fT; // underlying floating point type of type T
	using lT=typename cayley_hamilton<T>::lT; // longer bit representation of type T (if defined) or type T itself
	using cayley_hamilton<T>::n;
	using cayley_hamilton<T>::mmax;
	using cayley_hamilton<T>::nhl_max;
	using cayley_hamilton<T>::pl;
	using cayley_hamilton<T>::trpl;
	using cayley_hamilton<T>::crpl;
	using cayley_hamilton<T>::pal;
	using cayley_hamilton<T>::al;
	using cayley_hamilton<T>::tmat1;
	using cayley_hamilton<T>::tmat2;
	using cayley_hamilton<T>::set_n;
	using cayley_hamilton<T>::ch_mult;
	using cayley_hamilton<T>::set_to_zero;
	using cayley_hamilton<T>::set_to_identity;
	using cayley_hamilton<T>::set_to_identity_scaled;
	using cayley_hamilton<T>::matrix_copy;
	using cayley_hamilton<T>::matrix_copy_scaled;
	using cayley_hamilton<T>::matrix_add;
	using cayley_hamilton<T>::matrix_sub;
	using cayley_hamilton<T>::matrix_add_scaled;
	using cayley_hamilton<T>::matrix_mult_trace_nn;
	using cayley_hamilton<T>::matrix_mult_trace_na;
	using cayley_hamilton<T>::matrix_mult_trace_an;
	using cayley_hamilton<T>::matrix_mult_trace_aa;
	using cayley_hamilton<T>::matrix_mult_nn;
	using cayley_hamilton<T>::matrix_mult_na;
	using cayley_hamilton<T>::matrix_mult_an;
	using cayley_hamilton<T>::matrix_mult_scalar;
	using cayley_hamilton<T>::matrix_mult_scalar_add;
	using cayley_hamilton<T>::matrix_mult_scalar_sub;
	using cayley_hamilton<T>::matrix_mult_nn_add;
	using cayley_hamilton<T>::matrix_mult_na_add;
	using cayley_hamilton<T>::matrix_mult_an_add;
	using cayley_hamilton<T>::matrix_mult_nn_sub;
	using cayley_hamilton<T>::matrix_mult_na_sub;
	using cayley_hamilton<T>::matrix_mult_an_sub;
	using cayley_hamilton<T>::get_trace;
	using cayley_hamilton<T>::get_frobenius_norm;
	using cayley_hamilton<T>::get_one_norm;
	using cayley_hamilton<T>::get_l1_norm;
	using cayley_hamilton<T>::get_li_norm;
	using cayley_hamilton<T>::new_matrix;
	using cayley_hamilton<T>::delete_matrix;
	using cayley_hamilton<T>::new_matrix_array;
	using cayley_hamilton<T>::delete_matrix_array;
	using cayley_hamilton<T>::print_matrix;
	fT scalesquarelim;

	chexp() : cayley_hamilton<T>() {

	}

	chexp(int tn,fT sslim=1.0): cayley_hamilton<T>(tn),scalesquarelim(sslim) {

	}

	int operator()(T** ain,T** aout,int* onb=0) {
		// computes the matrix exponential of ain[][], using the Cayley-Hamilton recursion in combination with "scaling and squaring", 
		// and writes the result to aout[][]
		int niter=0;
		int nb=0;
		if(n>0) {
			int i,j,k;

			// determine nb for scaling factor 2^{-nb} of matrix ain[][] s.t. frobenius_norm(ain[][])*2^{-nb} \in (1.0,2.0) :
			fT anorm,sfac=1.0;
			get_frobenius_norm(ain,anorm);
			while(anorm*sfac>=scalesquarelim) {
				sfac*=0.5;
				++nb;
			}

			if(onb!=0) {
				*onb=nb;
			}

			//std::cout<<"nb="<<nb<<std::endl;

			matrix_copy_scaled(ain,sfac,pl[1]);

			trpl[0]=n;
			get_trace(pl[1],trpl[1]);
			for(i=2; i<n; ++i) {
				j=i/2;
				k=i%2;
				matrix_mult_nn(pl[j],pl[j+k],pl[i],trpl[i]);
			}
			j=n/2;
			k=n%2;
			matrix_mult_trace_nn(pl[j],pl[j+k],trpl[n]);

			// compute the characteristic polynomial crpl[] from the traced powers trpl[] :
			crpl[n]=1;
			for(j=1; j<=n; ++j) {
				crpl[n-j]=-trpl[j];
				for(i=1; i<j; ++i) {
					crpl[n-j]-=crpl[n-(j-i)]*trpl[i];
				}
				crpl[n-j]/=j;
			}

			// compute iteratively the n coefficients al[] so that the Cayley-Hamilton result
			// for the matrix power series is given by aout[][] = sum_i{al[i]*pl[i][][]} :

			// set initial values for the n entries in al[] and pal[] :
			fT wpf=1.0,twpf=1.0,wpfr=1.0,ttwpf; //leading coefficient of power series and of running sum used for convergence check
			for(i=0; i<n; ++i) {
				pal[i]=0;
				al[i]=wpf;
				wpf/=(fT)(i+1); //compute (i+1)-th power series coefficent from the i-th coefficient, using the rule defined by opf()
				twpf+=wpf;
			}
			pal[n-1]=1.0;

			// next we iteratively add higher order power series terms to al[] till al[] stops changing
			// more precisely: the iteration will terminate after nhl_max consecutive iterations leave twpf unchanged
			int nhl=0; // counts the number of consecutive non-changing iterations  
			T cho; // temporary variables for iterating
			fT s,rs=1.0; // used for normalizing the pal[]
			for(j=n; j<mmax; ++j) {
				s=0;
				cho=pal[n-1]*rs;
				for(i=n-1; i>0; --i) {
					pal[i]=pal[i-1]*rs-cho*crpl[i];
					s+=std::norm(pal[i]);
					al[i]+=wpf*pal[i];
				}
				pal[0]=-cho*crpl[0];
				s+=std::norm(pal[0]);
				al[0]+=wpf*pal[0];

				wpf/=(fT)(j+1); //compute (j+1)-th power series coefficent from j-th coefficient

				s=std::sqrt(s);
				if(s>1.0) {
					// if s is bigger than 1, normalize pal[] by a factor rs=1.0/s in next itaration, and multiply wpf by s to compensate
					wpf*=s;
					wpfr=1.0;
					rs=1.0/s;
				} else {
					wpfr=s;
					rs=1.0;
				}

				ttwpf=twpf;
				twpf+=wpf*wpfr;
				if(ttwpf==twpf) {
					++nhl;
					if(nhl>=nhl_max) {
						//terminate iteration
						break;
					}
				} else {
					nhl=0;
				}

			}
			niter=j;

			// take nb times the square of the result:
			for(k=0; k<nb; ++k) {
				ch_mult(crpl,al,al,al);
			}



			// form output matrix by summing the 0-th to (n-1)-th matrix powers pl[] with corresponding weights al[] 
			set_to_identity_scaled(al[0],aout);
			for(j=1; j<n; ++j) {
				matrix_add_scaled(pl[j],al[j],aout);
			}
		}
		return niter;
	}

	int operator()(T** ain,T** aout,T** daout) {
		// computes the matrix exponential of ain[][], using the Cayley-Hamilton recursion in combination with "scaling and squaring", 
		// and writes the result to aout[][]; computes also the derivative of aout[][] in the direction of
		// daout[][] and overwrite daout[][] with the result :
		int niter=0;
		int nb=0;
		if(n>0) {
			// use unused entries of pl[] as storage for computation of differentaial components:
			T** kmats=pl[0];
			T** kh=pl[n];

			int i,j,k,l;
			// determine nb for scaling factor 2^{-nb} of matrix ain[][] s.t. frobenius_norm(ain[][])*2^{-nb} \in (1.0,2.0) :
			fT anorm,sfac=1.0;
			get_frobenius_norm(ain,anorm);
			while(anorm*sfac>=scalesquarelim) {
				sfac*=0.5;
				++nb;
			}

			//std::cout<<"nb="<<nb<<std::endl;

			trpl[0]=n;

			matrix_copy_scaled(ain,sfac,pl[1]);

			get_trace(pl[1],trpl[1]);
			for(i=2; i<n; ++i) {
				j=i/2;
				k=i%2;
				matrix_mult_nn(pl[j],pl[j+k],pl[i],trpl[i]);
			}
			j=n/2;
			k=n%2;
			matrix_mult_trace_nn(pl[j],pl[j+k],trpl[n]);

			// compute the characteristic polynomial crpl[] from the traced powers trpl[] :
			crpl[n]=1;
			for(j=1; j<=n; ++j) {
				crpl[n-j]=-trpl[j];
				for(i=1; i<j; ++i) {
					crpl[n-j]-=crpl[n-(j-i)]*trpl[i];
				}
				crpl[n-j]/=j;
			}

			// compute iteratively the n coefficients al[] so that the Cayley-Hamilton result
			// for the matrix power series is given by aout[][] = sum_i{al[i]*pl[i][][]} :

			// set initial values for the n entries in al[] and pal[] and the nxn entries in kmats[][] :
			fT wpf=1.0,twpf=1.0,wpfr=1.0,ttwpf; //leading coefficient of power series and of running sum used for convergence check
			set_to_zero(kmats);
			for(i=0; i<n; ++i) {
				pal[i]=0;
				al[i]=wpf;
				wpf/=(fT)(i+1); //compute (i+1)-th power series coefficent from the i-th coefficient, using the rule defined by opf()

				k=i/2;
				for(j=i-k; j<=i; ++j) {
					kmats[i-j][j]=wpf;
				}
				twpf+=wpf;
			}
			pal[n-1]=1.0;
			// set initial values for the entries in kh[][] :
			set_to_zero(kh);
			k=(n-1)/2;
			for(i=n-1-k; i<n; ++i) {
				kh[n-1-i][i]=1.0;
			}


			// next we iteratively add higher order power series terms to al[] and kmats[][] till al[] stops changing
			// more precisely: the iteration will terminate after nhl_max consecutive iterations have not changed al[]	
			T cho; // temporary variables for iterating
			int nhl=0; // counts the number of consecutive non-changing iterations  
			fT s,rs=1.0; // used for normalizing the pal[]
			for(j=n; j<mmax; ++j) {
				s=0;
				cho=pal[n-1]*rs;
				for(i=n-1; i>0; --i) {
					pal[i]=pal[i-1]*rs-cho*crpl[i];
					s+=std::norm(pal[i]);
					al[i]+=wpf*pal[i];
				}
				pal[0]=-cho*crpl[0];
				s+=std::norm(pal[0]);
				al[0]+=wpf*pal[0];

				wpf/=(fT)(j+1); //compute (j+1)-th power series coefficent from j-th coefficient

				s=std::sqrt(s);
				if(s>1.0) {
					// if s is bigger than 1, normalize pal[] by a factor rs=1.0/s in next itaration, and multiply wpf by s to compensate
					wpf*=s;
					wpfr=1.0;
					rs=1.0/s;
				} else {
					wpfr=s;
					rs=1.0;
				}

				// add terms to kmats[][]:
				for(i=n-1; i>=0; --i) {
					cho=kh[i][n-1]*rs;
					for(k=n-1; k>i; --k) {
						kh[i][k]=kh[i][k-1]*rs-cho*crpl[k];
						kmats[i][k]+=wpf*kh[i][k];
					}
					if(i>0) {
						kh[i][i]=kh[i-1][i]*rs-cho*crpl[i];
					} else {
						kh[0][0]=pal[0]*rs-cho*crpl[0];
					}
					kmats[i][i]+=wpf*kh[i][i];
				}

				ttwpf=twpf;
				twpf+=wpf*wpfr;
				if(ttwpf==twpf) {
					++nhl;
					if(nhl>=nhl_max) {
						//terminate iteration
						break;
					}
				} else {
					nhl=0;
				}

			}
			niter=j;

			// take nb times the square:
			T* kal=trpl;
			for(k=0; k<nb; ++k) {
				cho=al[0];
				for(i=0; i<n; ++i) {
					kal[i]=pal[i]=al[i];
					kh[i][i]=(fT)0.5*kmats[i][i];
					kmats[i][i]=(fT)2.0*kh[0][i]*al[i];
					for(j=i+1; j<n; ++j) {
						kh[i][j]=(fT)0.5*kmats[i][j];
						kh[j][i]=(fT)0.5*kmats[i][j];
						kmats[i][j]=kh[0][i]*al[j]+kh[0][j]*al[i];
					}
					al[i]*=cho;
				}
				for(l=1; l<n; ++l) {
					cho=pal[n-1];
					for(i=n-1; i>0; --i) {
						pal[i]=pal[i-1]-cho*crpl[i];
						al[i]+=kal[l]*pal[i];
						for(j=i; j<n; ++j) {
							kmats[i][j]+=kh[i][l]*pal[j]+kh[l][j]*pal[i];
						}
					}
					pal[0]=-cho*crpl[0];
					al[0]+=kal[l]*pal[0];
					for(j=0; j<n; ++j) {
						kmats[0][j]+=kh[0][l]*pal[j]+kh[l][j]*pal[0];
					}
				}
			}


			// form output matrix by summing the 0-th to (n-1)-th matrix powers pl[] with corresponding weights al[] 
			set_to_identity_scaled(al[0],aout);
			for(j=1; j<n; ++j) {
				matrix_add_scaled(pl[j],al[j],aout);
			}


			matrix_copy(daout,tmat1);

			//i=0, j=0:
			set_to_identity_scaled(kmats[0][0],kh);
			//i=0, j>0:
			for(j=1; j<n; ++j) {
				matrix_add_scaled(pl[j],kmats[0][j],kh);
			}
			matrix_mult_nn(tmat1,kh,daout);

			//i>0:
			for(i=1; i<n; ++i) {
				//j=0:
				set_to_identity_scaled(kmats[0][i],kh);
				//j>0:
				for(j=1; j<i; ++j) {
					matrix_add_scaled(pl[j],kmats[j][i],kh);
				}
				for(j=i; j<n; ++j) {
					matrix_add_scaled(pl[j],kmats[i][j],kh);
				}
				matrix_mult_nn(tmat1,kh,tmat2);
				matrix_mult_nn_add(pl[i],tmat2,daout);
			}

		}
		return niter;
	}

	int operator()(T** ain,T** aout,T**** daout) {
		// computes the matrix exponential of ain[][], using the Cayley-Hamilton recursion in combination with "scaling and squaring", 
		// and writes the result to aout[][]; computes also the derivative of aout[][] with respect to each of 
		// the nxn components of ain[][] and write the result to daout[][][][] (the first two indices define the component
		// with respect to which the derivative is taken and the last two indices enumerate the components of the matrix-
		// valued derivative).
		int niter=0;
		int nb=0;
		if(n>0) {
			// use unused entries of pl[] as storage for computation of differentaial components:
			T** kmats=pl[0];
			T** kh=pl[n];

			int i,j,k,l;
			// determine nb for scaling factor 2^{-nb} of matrix ain[][] s.t. frobenius_norm(ain[][])*2^{-nb} \in (1.0,2.0) :
			fT anorm,sfac=1.0;
			get_frobenius_norm(ain,anorm);
			while(anorm*sfac>=scalesquarelim) {
				sfac*=0.5;
				++nb;
			}

			//std::cout<<"nb="<<nb<<std::endl;

			trpl[0]=n;

			matrix_copy_scaled(ain,sfac,pl[1]);

			get_trace(pl[1],trpl[1]);
			for(i=2; i<n; ++i) {
				j=i/2;
				k=i%2;
				matrix_mult_nn(pl[j],pl[j+k],pl[i],trpl[i]);
			}
			j=n/2;
			k=n%2;
			matrix_mult_trace_nn(pl[j],pl[j+k],trpl[n]);

			// compute the characteristic polynomial crpl[] from the traced powers trpl[] :
			crpl[n]=1;
			for(j=1; j<=n; ++j) {
				crpl[n-j]=-trpl[j];
				for(i=1; i<j; ++i) {
					crpl[n-j]-=crpl[n-(j-i)]*trpl[i];
				}
				crpl[n-j]/=j;
			}

			// compute iteratively the n coefficients al[] so that the Cayley-Hamilton result
			// for the matrix exponential is given by aout[][] = sum_i{al[i]*pl[i][][]} :

			// set initial values for the n entries in al[] and pal[] and the nxn entries in kmats[][] :
			fT wpf=1.0,twpf=1.0,wpfr=1.0,ttwpf; //leading coefficient of power series and of running sum used for convergence check
			set_to_zero(kmats);
			for(i=0; i<n; ++i) {
				pal[i]=0;
				al[i]=wpf;
				wpf/=(fT)(i+1); //compute (i+1)-th power series coefficent from the i-th coefficient, using the rule defined by opf()

				k=i/2;
				for(j=i-k; j<=i; ++j) {
					kmats[i-j][j]=wpf;
				}

				twpf+=wpf;
			}
			pal[n-1]=1.0;
			// set initial values for the entries in kh[][] :
			set_to_zero(kh);
			k=(n-1)/2;
			for(i=n-1-k; i<n; ++i) {
				kh[n-1-i][i]=1.0;
			}


			// next we iteratively add higher order power series terms to al[] and kmats[][] till al[] stops changing
			// more precisely: the iteration will terminate after nhl_max consecutive iterations have not changed al[]	
			T cho; // temporary variables for iterating
			int nhl=0; // counts the number of consecutive non-changing iterations  
			fT s,rs=1.0; // used for normalizing the pal[]
			for(j=n; j<mmax; ++j) {
				s=0;
				cho=pal[n-1]*rs;
				for(i=n-1; i>0; --i) {
					pal[i]=pal[i-1]*rs-cho*crpl[i];
					s+=std::norm(pal[i]);
					al[i]+=wpf*pal[i];
				}
				pal[0]=-cho*crpl[0];
				s+=std::norm(pal[0]);
				al[0]+=wpf*pal[0];

				wpf/=(fT)(j+1); //compute (j+1)-th power series coefficent from j-th coefficient

				s=std::sqrt(s);
				if(s>1.0) {
					// if s is bigger than 1, normalize pal[] by a factor rs=1.0/s in next itaration, and multiply wpf by s to compensate
					wpf*=s;
					wpfr=1.0;
					rs=1.0/s;
				} else {
					wpfr=s;
					rs=1.0;
				}

				// add terms to kmats[][]:
				for(i=n-1; i>=0; --i) {
					cho=kh[i][n-1]*rs;
					for(k=n-1; k>i; --k) {
						kh[i][k]=kh[i][k-1]*rs-cho*crpl[k];
						kmats[i][k]+=wpf*kh[i][k];
					}
					if(i>0) {
						kh[i][i]=kh[i-1][i]*rs-cho*crpl[i];
					} else {
						kh[0][0]=pal[0]*rs-cho*crpl[0];
					}
					kmats[i][i]+=wpf*kh[i][i];
				}

				ttwpf=twpf;
				twpf+=wpf*wpfr;
				if(ttwpf==twpf) {
					++nhl;
					if(nhl>=nhl_max) {
						//terminate iteration
						break;
					}
				} else {
					nhl=0;
				}

			}
			niter=j;

			// compute coefficients for nb-times squared matrix function and corresponding differentials:
			T* kal=trpl;
			for(k=0; k<nb; ++k) {
				cho=al[0];
				for(i=0; i<n; ++i) {
					kal[i]=pal[i]=al[i];
					kh[i][i]=(fT)0.5*kmats[i][i];
					kmats[i][i]=(fT)2.0*kh[0][i]*al[i];
					for(j=i+1; j<n; ++j) {
						kh[j][i]=kh[i][j]=(fT)0.5*kmats[i][j]; // define symmetric kh[][] to avoid case-distinctions in computations below
						kmats[i][j]=kh[0][i]*al[j]+kh[0][j]*al[i];
					}
					al[i]*=cho;
				}
				for(l=1; l<n; ++l) {
					cho=pal[n-1];
					for(i=n-1; i>0; --i) {
						pal[i]=pal[i-1]-cho*crpl[i];
						al[i]+=kal[l]*pal[i];
						for(j=i; j<n; ++j) {
							kmats[i][j]+=kh[i][l]*pal[j]+kh[l][j]*pal[i];
						}
					}
					pal[0]=-cho*crpl[0];
					al[0]+=kal[l]*pal[0];
					for(j=0; j<n; ++j) {
						kmats[0][j]+=kh[0][l]*pal[j]+kh[l][j]*pal[0];
					}
				}
			}

			// form output matrix by summing the 0-th to (n-1)-th matrix powers pl[] with corresponding weights al[] 
			set_to_identity_scaled(al[0],aout);
			for(j=1; j<n; ++j) {
				matrix_add_scaled(pl[j],al[j],aout);
			}

			// form output for matrix of differentials daout[ic1][ic2][k][l] = daout[k][l]/dain[ic2][ic1] :
			// (note: in principle one could use symmetry domat[ic1][ic2].e(k, l) = domat[l][k].e(ic2, ic1),
			//  which would amount to  (i, j) <--> (j, i) with i = ic1 + n * ic2;  j = l + n * k; )

			//i=0, j=0:
			set_to_identity_scaled(kmats[0][0],kh);
			//i=0, j>0:
			for(j=1; j<n; ++j) {
				matrix_add_scaled(pl[j],kmats[0][j],kh);
			}
			int ic1,ic2;
			T** tdaout;
			for(ic1=0; ic1<n; ++ic1) {
				for(ic2=0; ic2<n; ++ic2) {
					tdaout=daout[ic1][ic2];
					set_to_zero(tdaout);
					for(k=0; k<n; ++k) {
						tdaout[k][ic1]+=kh[k][ic2];
					}
				}
			}

			//i>0:
			for(i=1; i<n; ++i) {
				//j=0:
				set_to_identity_scaled(kmats[0][i],kh);
				//j>0:
				for(j=1; j<i; ++j) {
					matrix_add_scaled(pl[j],kmats[j][i],kh);
				}
				for(j=i; j<n; ++j) {
					matrix_add_scaled(pl[j],kmats[i][j],kh);
				}
				for(ic1=0; ic1<n; ++ic1) {
					for(ic2=0; ic2<n; ++ic2) {
						tdaout=daout[ic1][ic2];
						for(k=0; k<n; ++k) {
							for(l=0; l<n; ++l) {
								tdaout[k][l]+=pl[i][ic1][l]*kh[k][ic2];
							}
						}

					}
				}

			}

		}
		return niter;
	}

	int get_r_k(T** ain,T* rout,T** kout) {
		// computes the matrix exponential of ain[][], using the Cayley-Hamilton recursion in combination with "scaling and squaring", 
		// and writes the Cayley-Hamilton coefficients to rout[]; computes also the decomposition matrix of the derivative
		// of rout[] with respect to the components of ain[][] and writes the result to kout[][] :
		int niter=0;
		int nb=0;
		if(n>0) {
			T* rl=rout;
			// use unused entries of pl[] as storage for computation of differentaial components:
			T** kmats=kout;
			T** kh=pl[n];

			int i,j,k,l;
			// determine nb for scaling factor 2^{-nb} of matrix ain[][] s.t. frobenius_norm(ain[][])*2^{-nb} \in (1.0,2.0) :
			fT anorm,sfac=1.0;
			get_frobenius_norm(ain,anorm);
			while(anorm*sfac>=scalesquarelim) {
				sfac*=0.5;
				++nb;
			}

			//std::cout<<"nb="<<nb<<std::endl;

			trpl[0]=n;

			matrix_copy_scaled(ain,sfac,pl[1]);

			get_trace(pl[1],trpl[1]);
			for(i=2; i<n; ++i) {
				j=i/2;
				k=i%2;
				matrix_mult_nn(pl[j],pl[j+k],pl[i],trpl[i]);
			}
			j=n/2;
			k=n%2;
			matrix_mult_trace_nn(pl[j],pl[j+k],trpl[n]);

			// compute the characteristic polynomial crpl[] from the traced powers trpl[] :
			crpl[n]=1;
			for(j=1; j<=n; ++j) {
				crpl[n-j]=-trpl[j];
				for(i=1; i<j; ++i) {
					crpl[n-j]-=crpl[n-(j-i)]*trpl[i];
				}
				crpl[n-j]/=j;
			}

			// compute iteratively the n coefficients al[] so that the Cayley-Hamilton result
			// for the matrix exponential would be given by aout[][] = sum_i{al[i]*pl[i][][]} :

			// set initial values for the n entries in al[] and pal[] and the nxn entries in kmats[][] :
			fT wpf=1.0,twpf=1.0,wpfr=1.0,ttwpf; //leading coefficient of power series and of running sum used for convergence check
			set_to_zero(kmats);
			for(i=0; i<n; ++i) {
				pal[i]=0;
				rl[i]=wpf;
				wpf/=(fT)(i+1); //compute (i+1)-th power series coefficent from the i-th coefficient, using the rule defined by opf()

				k=i/2;
				for(j=i-k; j<=i; ++j) {
					kmats[i-j][j]=wpf;
				}

				twpf+=wpf;
			}
			pal[n-1]=1.0;
			// set initial values for the entries in kh[][] :
			set_to_zero(kh);
			k=(n-1)/2;
			for(i=n-1-k; i<n; ++i) {
				kh[n-1-i][i]=1.0;
			}


			// next we iteratively add higher order power series terms to al[] and kmats[][] till al[] stops changing
			// more precisely: the iteration will terminate after nhl_max consecutive iterations have not changed al[]	
			T cho; // temporary variables for iterating
			int nhl=0; // counts the number of consecutive non-changing iterations  
			fT s,rs=1.0; // used for normalizing the pal[]
			for(j=n; j<mmax; ++j) {
				s=0;
				cho=pal[n-1]*rs;
				for(i=n-1; i>0; --i) {
					pal[i]=pal[i-1]*rs-cho*crpl[i];
					s+=std::norm(pal[i]);
					rl[i]+=wpf*pal[i];
				}
				pal[0]=-cho*crpl[0];
				s+=std::norm(pal[0]);
				rl[0]+=wpf*pal[0];

				wpf/=(fT)(j+1); //compute (j+1)-th power series coefficent from j-th coefficient

				s=std::sqrt(s);
				if(s>1.0) {
					// if s is bigger than 1, normalize pal[] by a factor rs=1.0/s in next itaration, and multiply wpf by s to compensate
					wpf*=s;
					wpfr=1.0;
					rs=1.0/s;
				} else {
					wpfr=s;
					rs=1.0;
				}

				// add terms to kmats[][]:
				for(i=n-1; i>=0; --i) {
					cho=kh[i][n-1]*rs;
					for(k=n-1; k>i; --k) {
						kh[i][k]=kh[i][k-1]*rs-cho*crpl[k];
						kmats[i][k]+=wpf*kh[i][k];
					}
					if(i>0) {
						kh[i][i]=kh[i-1][i]*rs-cho*crpl[i];
					} else {
						kh[0][0]=pal[0]*rs-cho*crpl[0];
					}
					kmats[i][i]+=wpf*kh[i][i];
				}

				ttwpf=twpf;
				twpf+=wpf*wpfr;
				if(ttwpf==twpf) {
					++nhl;
					if(nhl>=nhl_max) {
						//terminate iteration
						break;
					}
				} else {
					nhl=0;
				}

			}
			niter=j;

			// take nb times the square of the result matrix:
			T* kal=trpl;
			for(k=0; k<nb; ++k) {
				cho=rl[0];
				for(i=0; i<n; ++i) {
					kal[i]=pal[i]=rl[i];
					kh[i][i]=(fT)0.5*kmats[i][i];
					kmats[i][i]=(fT)2.0*kh[0][i]*rl[i];
					for(j=i+1; j<n; ++j) {
						kh[i][j]=(fT)0.5*kmats[i][j];
						kh[j][i]=(fT)0.5*kmats[i][j];
						kmats[i][j]=kh[0][i]*rl[j]+kh[0][j]*rl[i];
					}
					rl[i]*=cho;
				}
				for(l=1; l<n; ++l) {
					cho=pal[n-1];
					for(i=n-1; i>0; --i) {
						pal[i]=pal[i-1]-cho*crpl[i];
						rl[i]+=kal[l]*pal[i];
						for(j=i; j<n; ++j) {
							kmats[i][j]+=kh[i][l]*pal[j]+kh[l][j]*pal[i];
						}
					}
					pal[0]=-cho*crpl[0];
					rl[0]+=kal[l]*pal[0];
					for(j=0; j<n; ++j) {
						kmats[0][j]+=kh[0][l]*pal[j]+kh[l][j]*pal[0];
					}
				}
			}
			if(nb>0) {
				// input matrix has been scaled :
				// 
				// rescale elements of rl[] to get coefficients w.r.t. to powers of originial input matrix ain[][] :
				s=sfac;
				for(i=1; i<n; ++i) {
					rl[i]*=s;
					s*=sfac;
				}
				// rescale elements of kmats[][] to get coefficients w.r.t. to powers of originial input matrix ain[][] :
				s=1.0;
				for(i=0; i<n; ++i) {
					rs=s;
					for(j=i; j<n; ++j) {
						kmats[i][j]*=rs;
						rs*=sfac;
					}
					s*=sfac*sfac; // (starting each row from diagonal)
				}
			}

		}
		return niter;
	}

	int get_r_dr(T** ain,T* rout,T*** drout) {
		// computes the matrix exponential of ain[][], using the Cayley-Hamilton recursion in combination with "scaling and squaring", 
		// and writes the Cayley-Hamilton coefficients to rout[]; computes also the derivatives of rout[]
		// with respect to the components of ain[][] and writes the resulting matrices to drout[][][] :
		int niter=0;
		int nb=0;
		if(n>0) {
			T* rl=rout;
			// use unused entries of pl[] as storage for computation of differentaial components:
			T** kmats=pl[0];
			T** kh=pl[n];

			int i,j,k,l;
			// determine nb for scaling factor 2^{-nb} of matrix ain[][] s.t. frobenius_norm(ain[][])*2^{-nb} \in (1.0,2.0) :
			fT anorm,sfac=1.0;
			get_frobenius_norm(ain,anorm);
			while(anorm*sfac>=scalesquarelim) {
				sfac*=0.5;
				++nb;
			}

			//std::cout<<"nb="<<nb<<std::endl;

			trpl[0]=n;

			matrix_copy_scaled(ain,sfac,pl[1]);

			get_trace(pl[1],trpl[1]);
			for(i=2; i<n; ++i) {
				j=i/2;
				k=i%2;
				matrix_mult_nn(pl[j],pl[j+k],pl[i],trpl[i]);
			}
			j=n/2;
			k=n%2;
			matrix_mult_trace_nn(pl[j],pl[j+k],trpl[n]);

			// compute the characteristic polynomial crpl[] from the traced powers trpl[] :
			crpl[n]=1;
			for(j=1; j<=n; ++j) {
				crpl[n-j]=-trpl[j];
				for(i=1; i<j; ++i) {
					crpl[n-j]-=crpl[n-(j-i)]*trpl[i];
				}
				crpl[n-j]/=j;
			}

			// compute iteratively the n coefficients al[] so that the Cayley-Hamilton result
			// for the matrix exponential would be given by aout[][] = sum_i{al[i]*pl[i][][]} :

			// set initial values for the n entries in al[] and pal[] and the nxn entries in kmats[][] :
			fT wpf=1.0,twpf=1.0,wpfr=1.0,ttwpf; //leading coefficient of power series and of running sum used for convergence check
			set_to_zero(kmats);
			for(i=0; i<n; ++i) {
				pal[i]=0;
				rl[i]=wpf;
				wpf/=(fT)(i+1); //compute (i+1)-th power series coefficent from the i-th coefficient, using the rule defined by opf()

				k=i/2;
				for(j=i-k; j<=i; ++j) {
					kmats[i-j][j]=wpf;
				}

				twpf+=wpf;
			}
			pal[n-1]=1.0;
			// set initial values for the entries in kh[][] :
			set_to_zero(kh);
			k=(n-1)/2;
			for(i=n-1-k; i<n; ++i) {
				kh[n-1-i][i]=1.0;
			}


			// next we iteratively add higher order power series terms to al[] and kmats[][] till al[] stops changing
			// more precisely: the iteration will terminate after nhl_max consecutive iterations have not changed al[]	
			T cho; // temporary variables for iterating
			int nhl=0; // counts the number of consecutive non-changing iterations  
			fT s,rs=1.0; // used for normalizing the pal[]
			for(j=n; j<mmax; ++j) {
				s=0;
				cho=pal[n-1]*rs;
				for(i=n-1; i>0; --i) {
					pal[i]=pal[i-1]*rs-cho*crpl[i];
					s+=std::norm(pal[i]);
					rl[i]+=wpf*pal[i];
				}
				pal[0]=-cho*crpl[0];
				s+=std::norm(pal[0]);
				rl[0]+=wpf*pal[0];

				wpf/=(fT)(j+1); //compute (j+1)-th power series coefficent from j-th coefficient

				s=std::sqrt(s);
				if(s>1.0) {
					// if s is bigger than 1, normalize pal[] by a factor rs=1.0/s in next itaration, and multiply wpf by s to compensate
					wpf*=s;
					wpfr=1.0;
					rs=1.0/s;
				} else {
					wpfr=s;
					rs=1.0;
				}

				// add terms to kmats[][]:
				for(i=n-1; i>=0; --i) {
					cho=kh[i][n-1]*rs;
					for(k=n-1; k>i; --k) {
						kh[i][k]=kh[i][k-1]*rs-cho*crpl[k];
						kmats[i][k]+=wpf*kh[i][k];
					}
					if(i>0) {
						kh[i][i]=kh[i-1][i]*rs-cho*crpl[i];
					} else {
						kh[0][0]=pal[0]*rs-cho*crpl[0];
					}
					kmats[i][i]+=wpf*kh[i][i];
				}

				ttwpf=twpf;
				twpf+=wpf*wpfr;
				if(ttwpf==twpf) {
					++nhl;
					if(nhl>=nhl_max) {
						//terminate iteration
						break;
					}
				} else {
					nhl=0;
				}

			}
			niter=j;

			// take nb times the square of the result matrix:
			T* kal=trpl;
			for(k=0; k<nb; ++k) {
				cho=rl[0];
				for(i=0; i<n; ++i) {
					kal[i]=pal[i]=rl[i];
					kh[i][i]=(fT)0.5*kmats[i][i];
					kmats[i][i]=(fT)2.0*kh[0][i]*rl[i];
					for(j=i+1; j<n; ++j) {
						kh[i][j]=(fT)0.5*kmats[i][j];
						kh[j][i]=(fT)0.5*kmats[i][j];
						kmats[i][j]=kh[0][i]*rl[j]+kh[0][j]*rl[i];
					}
					rl[i]*=cho;
				}
				for(l=1; l<n; ++l) {
					cho=pal[n-1];
					for(i=n-1; i>0; --i) {
						pal[i]=pal[i-1]-cho*crpl[i];
						rl[i]+=kal[l]*pal[i];
						for(j=i; j<n; ++j) {
							kmats[i][j]+=kh[i][l]*pal[j]+kh[l][j]*pal[i];
						}
					}
					pal[0]=-cho*crpl[0];
					rl[0]+=kal[l]*pal[0];
					for(j=0; j<n; ++j) {
						kmats[0][j]+=kh[0][l]*pal[j]+kh[l][j]*pal[0];
					}
				}
			}

			// compute drout[][][] from the kmats[][] matrix :

			// compute \tilde{R}-matrix:
			for(i=0; i<n; ++i) {
				tmat1[i][0]=kmats[i][n-1];
			}
			for(j=1; j<n; ++j) {
				cho=tmat1[n-1][j-1];
				for(i=n-1; i>0; --i) {
					tmat1[i][j]=tmat1[i-1][j-1]-cho*crpl[i];
				}
				tmat1[0][j]=-cho*crpl[0];
			}

			// compute R-matrix:
			for(j=0; j<n; ++j) {
				for(i=0; i<n; ++i) {
					tmat2[j][i]=tmat1[j][n-i-1];
					for(k=0; k<n-i-1; ++k) {
						tmat2[j][i]+=crpl[1+i+k]*tmat1[j][k];
					}
				}
			}
			// compute drout[]-matrices:
			for(k=0; k<n; ++k) {
				set_to_identity_scaled(tmat2[k][0],drout[k]);
				for(i=1; i<n; ++i) {
					matrix_mult_scalar_add(pl[i],tmat2[k][i],drout[k]);
				}
			}

			if(nb>0) {
				// input matrix has been scaled :

				// rescale elements of rl[] and drout[][][] to get coefficients w.r.t. to powers of originial input matrix ain[][] :
				s=sfac;
				for(i=1; i<n; ++i) {
					rl[i]*=s;
					matrix_mult_scalar(drout[i],s,drout[i]);
					s*=sfac;
				}
			}

		}
		return niter;
	}

};

template<class T>
class nvexp: public cayley_hamilton<T> {
// template class providing an implementation of matrix exponentiation using naive Taylor series method
// with scaling and squaring (cf. arXiv:2404.07704).
public:
	using fT=typename cayley_hamilton<T>::fT; // underlying floating point type of type T
	using lT=typename cayley_hamilton<T>::lT; // longer bit representation of type T (if defined) or type T itself
	using cayley_hamilton<T>::n;
	using cayley_hamilton<T>::mmax;
	using cayley_hamilton<T>::nhl_max;
	using cayley_hamilton<T>::pl;
	using cayley_hamilton<T>::trpl;
	using cayley_hamilton<T>::crpl;
	using cayley_hamilton<T>::pal;
	using cayley_hamilton<T>::al;
	using cayley_hamilton<T>::tmat1;
	using cayley_hamilton<T>::tmat2;
	using cayley_hamilton<T>::set_n;
	using cayley_hamilton<T>::ch_mult;
	using cayley_hamilton<T>::set_to_zero;
	using cayley_hamilton<T>::set_to_identity;
	using cayley_hamilton<T>::set_to_identity_scaled;
	using cayley_hamilton<T>::matrix_copy;
	using cayley_hamilton<T>::matrix_copy_scaled;
	using cayley_hamilton<T>::matrix_add;
	using cayley_hamilton<T>::matrix_sub;
	using cayley_hamilton<T>::matrix_add_scaled;
	using cayley_hamilton<T>::matrix_mult_trace_nn;
	using cayley_hamilton<T>::matrix_mult_trace_na;
	using cayley_hamilton<T>::matrix_mult_trace_an;
	using cayley_hamilton<T>::matrix_mult_trace_aa;
	using cayley_hamilton<T>::matrix_mult_nn;
	using cayley_hamilton<T>::matrix_mult_na;
	using cayley_hamilton<T>::matrix_mult_an;
	using cayley_hamilton<T>::matrix_mult_scalar;
	using cayley_hamilton<T>::matrix_mult_scalar_add;
	using cayley_hamilton<T>::matrix_mult_scalar_sub;
	using cayley_hamilton<T>::matrix_mult_nn_add;
	using cayley_hamilton<T>::matrix_mult_na_add;
	using cayley_hamilton<T>::matrix_mult_an_add;
	using cayley_hamilton<T>::matrix_mult_nn_sub;
	using cayley_hamilton<T>::matrix_mult_na_sub;
	using cayley_hamilton<T>::matrix_mult_an_sub;
	using cayley_hamilton<T>::get_trace;
	using cayley_hamilton<T>::get_frobenius_norm;
	using cayley_hamilton<T>::get_one_norm;
	using cayley_hamilton<T>::get_l1_norm;
	using cayley_hamilton<T>::get_li_norm;
	using cayley_hamilton<T>::new_matrix;
	using cayley_hamilton<T>::delete_matrix;
	using cayley_hamilton<T>::new_matrix_array;
	using cayley_hamilton<T>::delete_matrix_array;
	using cayley_hamilton<T>::print_matrix;
	fT scalesquarelim;

	nvexp(): cayley_hamilton<T>() {

	}

	nvexp(int tn, fT sslim=1.0): cayley_hamilton<T>(tn),scalesquarelim(sslim) {

	}

	int operator()(T** ain,T** aout,int* onb=0) {
		// computes the matrix exponential of ain[][], using naive Taylor series method with "scaling and squaring", 
		// and writes the result to aout[][]
		int niter=0;
		int nb=0;
		if(n>0) {
			int j,k;

			// determine nb for scaling factor 2^{-nb} of matrix ain[][] s.t. frobenius_norm(ain[][])*2^{-nb} \in (1.0,2.0) :
			fT anorm,sfac=1.0;
			get_frobenius_norm(ain,anorm);
			while(anorm*sfac>=scalesquarelim) {
				sfac*=0.5;
				++nb;
			}

			if(onb!=0) {
				*onb=nb;
			}

			//std::cout<<"nb="<<nb<<std::endl;

			matrix_copy_scaled(ain,sfac,tmat1);
			matrix_copy(tmat1,pl[0]);
			set_to_identity(aout);
			matrix_add(pl[0],aout);
			int nhl=0; // counts the number of consecutive non-changing iterations  
			fT ttwpf,twpf;
			get_frobenius_norm(aout,twpf);
			for(j=2; j<mmax; ++j) {
				matrix_mult_nn(tmat1,pl[0],pl[1]);
				matrix_mult_scalar(pl[1],(fT)1.0/(fT)j,pl[0]);
				matrix_add(pl[0],aout);
				ttwpf=twpf;
				get_frobenius_norm(aout,twpf);
				if(ttwpf==twpf) {
					++nhl;
					if(nhl>=nhl_max) {
						//terminate iteration
						break;
					}
				} else {
					nhl=0;
				}
			}
			niter=j;

			// take nb times the square of the result:
			T* kal=trpl;
			for(k=0; k<nb; ++k) {
				matrix_mult_nn(aout,aout,pl[0]);
				matrix_copy(pl[0],aout);
			}

		}
		return niter;
	}

};