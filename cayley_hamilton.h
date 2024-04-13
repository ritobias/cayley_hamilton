/*
	header-only C++ implementation of the iterative Cayley-Hamilton method 
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

typedef double ftype; // floating point representation to be used
typedef long double lftype; // long version of ftype

typedef std::complex<ftype> ctype; // complex type corresponding to ftype
typedef std::complex<lftype> lctype; // complex type corresponding to lftype

typedef void (*opf)(ftype*,ftype*,ftype*);

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

	cayley_hamilton(): n(0),a(0),pl(0),trpl(0),crpl(0),dpl(0),dtrpl(0),dcrpl(0),pal(0),dpal(0),al(0),dal(0),mmax(0),nhl_max(0),tmat(0),wps(0) {
		//default constructor; will require a call to set_n() to set the size of the square matrices on which the class will operate

		opf=[](fT pref,int i) { return pref/(fT)i; }; // returns the i-th coeff. of the power seris, computed from the (i-1)-th coeff. "pref" and "i"
													        // default rule generates the coeffs. for the powerseries of exp(), i.e. 1/i!
	}

	cayley_hamilton(int tn,fT(*topf)(fT,int)=0) {
		// constrator with matrix size n=tn as argument and optional argument for function pointer that defines power series coeffients. 
		if(topf!=0) {
			opf=topf;
		}
		if(tn>0) {
			n=tn;
			new_matrix(a);
			new_matrix(tmat);
			new_matrix_array(pl,n+1);
			dpl=new T***[n*n];
			for(int i=0; i<n*n; ++i) {
				new_matrix_array(dpl[i],n+1);
			}
			trpl=new T[n+1]();
			dtrpl=new T*[n*n];
			for(int i=0; i<n*n; ++i) {
				dtrpl[i]=new T[n+1]();
			}
			crpl=new T[n+1]();
			dcrpl=new T*[n*n];
			for(int i=0; i<n*n; ++i) {
				dcrpl[i]=new T[n+1]();
			}
			pal=new T[n]();
			dpal=new T*[n*n];
			for(int i=0; i<n*n; ++i) {
				dpal[i]=new T[n]();
			}
			al=new T[n]();
			dal=new T*[n*n];
			for(int i=0; i<n*n; ++i) {
				dal[i]=new T[n]();
			}
			wps=new fT[n]();
			opf=[](fT pref,int i) { return pref/(fT)i; };
			mmax=100*n;
			nhl_max=3;
		} else {
			n=0;
			a=0;
			tmat=0;
			pl=0;
			dpl=0;
			trpl=0;
			dtrpl=0;
			crpl=0;
			dcrpl=0;
			pal=0;
			dpal=0;
			al=0;
			dal=0;
			wps=0;
			opf=[](fT pref,int i) { return pref/(fT)i; };
			mmax=100;
			nhl_max=3;
		}
	}

	~cayley_hamilton() {
		if(a!=0) {
			delete_matrix(a);
		}
		if(tmat!=0) {
			delete_matrix(tmat);
		}
		if(pl!=0) {
			delete_matrix_array(pl,n+1);
		}
		if(dpl!=0) {
			for(int i=0; i<n*n; ++i) {
				delete_matrix_array(dpl[i],n+1);
			}
			delete[] dpl;
			dpl=0;
		}
		if(trpl!=0) {
			delete[] trpl;
			trpl=0;
		}
		if(dtrpl!=0) {
			for(int i=0; i<n*n; ++i) {
				delete[] dtrpl[i];
			}
			delete[] dtrpl;
			dtrpl=0;
		}
		if(crpl!=0) {
			delete[] crpl;
			crpl=0;
		}
		if(dcrpl!=0) {
			for(int i=0; i<n*n; ++i) {
				delete[] dcrpl[i];
			}
			delete[] dcrpl;
			dcrpl=0;
		}
		if(pal!=0) {
			delete[] pal;
			pal=0;
		}
		if(dpal!=0) {
			for(int i=0; i<n*n; ++i) {
				delete[] dpal[i];
			}
			delete[] dpal;
			dpal=0;
		}
		if(al!=0) {
			delete[] al;
			al=0;
		}
		if(dal!=0) {
			for(int i=0; i<n*n; ++i) {
				delete[] dal[i];
			}
			delete[] dal;
			dal=0;
		}
		if(wps!=0) {
			delete[] wps;
			wps=0;
		}
		n=0;
	}

	void set_n(int tn) {
		//if an instance has been created with the trivial/empty constructor, 
		//this function has to be run afterwards to initialize the instance properly
		//and to prepare it to operate on matrices of size nxn with n=tn
		if(tn>0) {
			if(tn!=n) {
				if(a!=0) {
					delete_matrix(a);
				}
				if(tmat!=0) {
					delete_matrix(tmat);
				}
				if(pl!=0) {
					delete_matrix_array(pl,n+1);
				}
				if(dpl!=0) {
					for(int i=0; i<n*n; ++i) {
						delete_matrix_array(dpl[i],n+1);
					}
					delete[] dpl;
					dpl=0;
				}
				if(trpl!=0) {
					delete[] trpl;
					trpl=0;
				}
				if(dtrpl!=0) {
					for(int i=0; i<n*n; ++i) {
						delete[] dtrpl[i];
					}
					delete[] dtrpl;
					dtrpl=0;
				}
				if(crpl!=0) {
					delete[] crpl;
					crpl=0;
				}
				if(dcrpl!=0) {
					for(int i=0; i<n*n; ++i) {
						delete[] dcrpl[i];
					}
					delete[] dcrpl;
					dcrpl=0;
				}
				if(pal!=0) {
					delete[] pal;
					pal=0;
				}
				if(dpal!=0) {
					for(int i=0; i<n*n; ++i) {
						delete[] dpal[i];
					}
					delete[] dpal;
					dpal=0;
				}
				if(al!=0) {
					delete[] al;
					al=0;
				}
				if(dal!=0) {
					for(int i=0; i<n*n; ++i) {
						delete[] dal[i];
					}
					delete[] dal;
					dal=0;
				}
				if(wps!=0) {
					delete[] wps;
					wps=0;
				}

				n=tn;

				new_matrix(a);

				new_matrix(tmat);

				new_matrix_array(pl,n+1);

				dpl=new T***[n*n];
				for(int i=0; i<n*n; ++i) {
					new_matrix_array(dpl[i],n+1);
				}

				trpl=new T[n+1];

				dtrpl=new T*[n*n];
				for(int i=0; i<n*n; ++i) {
					dtrpl[i]=new T[n+1]();
				}

				crpl=new T[n+1];

				dcrpl=new T*[n*n];
				for(int i=0; i<n*n; ++i) {
					dcrpl[i]=new T[n+1]();
				}

				pal=new T[n];

				dpal=new T*[n*n];
				for(int i=0; i<n*n; ++i) {
					dpal[i]=new T[n+1]();
				}

				al=new T[n];

				dal=new T*[n*n];
				for(int i=0; i<n*n; ++i) {
					dal[i]=new T[n+1]();
				}
				
				wps=new fT[n]();

				mmax=100*n;
			}
		} else {
			if(a!=0) {
				delete_matrix(a);
			}
			if(tmat!=0) {
				delete_matrix(tmat);
			}
			if(pl!=0) {
				delete_matrix_array(pl,n+1);
			}
			if(dpl!=0) {
				for(int i=0; i<n*n; ++i) {
					delete_matrix_array(dpl[i],n+1);
				}
				delete[] dpl;
				dpl=0;
			}
			if(trpl!=0) {
				delete[] trpl;
				trpl=0;
			}
			if(dtrpl!=0) {
				for(int i=0; i<n*n; ++i) {
					delete[] dtrpl[i];
				}
				delete[] dtrpl;
				dtrpl=0;
			}
			if(crpl!=0) {
				delete[] crpl;
				crpl=0;
			}
			if(dcrpl!=0) {
				for(int i=0; i<n*n; ++i) {
					delete[] dcrpl[i];
				}
				delete[] dcrpl;
				dcrpl=0;
			}
			if(pal!=0) {
				delete[] pal;
				pal=0;
			}
			if(dpal!=0) {
				for(int i=0; i<n*n; ++i) {
					delete[] dpal[i];
				}
				delete[] dpal;
				dpal=0;
			}
			if(al!=0) {
				delete[] al;
				al=0;
			}
			if(dal!=0) {
				for(int i=0; i<n*n; ++i) {
					delete[] dal[i];
				}
				delete[] dal;
				dal=0;
			}
			if(wps!=0) {
				delete[] wps;
				wps=0;
			}
			mmax=0;

			n=0;
		}
	}

	void set_opf(fT(*topf)(fT,int)) {
		//set opf to point to a user-defined function or lambda for generating the power series coefficients
		opf=topf;
	}

	void operator()(T** ain,T** aout,int rescale=0) {
		// applies the power series defined by opf() to the matrix ain[][], using the Cayley-Hamilton recursion, 
		// and writes the result to aout[][]
		// if rescale is set to 1, then the computation will be performed with rescaled input matrix, which is
		// useful for matrices of large norm; if rescale isset to 0, no matrix rescaling will be performed. 

		int i,j,k;
		fT sfac=1.0; //scaling factor
		// compute the 0-th to n-th matrix powers of ta[][] :
		//   the i-th matrix power of ta[][] is stored in pl[i][][]
		//   the trace of the i-th matrix power of ta[][] is stored in trpl[i]
		set_to_identity(pl[0]);
		trpl[0]=n;
		if(rescale==1) {
			// if matrix rescaling is used, set sfac to be the magnitude of the largest element of ain[][]
			// and initiate the computation of the matrix powers from ain[][]/sfac instead of ain[][]:
			// (note that since we compute also the matrix powers pl[] from the rescaled matrix,
			// the Cayley-Hamilton coefficient a_{k,j}, with k=0,1,...,\infty; j=0,...,n-1 
			// will need to be rescaled only by a factor  sfac^{k}, instead of sfac^{k-j})
			sfac=max_matrix_el(ain);
			matrix_copy_scaled(ain,1.0/sfac,pl[1]);
		} else {
			matrix_copy(ain,pl[1]);
		}
		get_trace(pl[1],trpl[1]);
		for(i=2; i<=n; ++i) {
			j=i/2;
			k=i%2;
			matrix_mult_nn(pl[j],pl[j+k],pl[i],trpl[i]);
		}

		// compute the characteristic polynomial crpl[] from the traced powers trpl[] :
		crpl[n]=1;
		for(j=1; j<=n; ++j) {
			crpl[n-j]=0;
			for(i=1; i<=j; ++i) {
				crpl[n-j]-=crpl[n-(j-i)]*trpl[i];
			}
			crpl[n-j]/=j;
		}

		// compute iteratively the n coefficients al[] so that the Cayley-Hamilton result
		// for the matrix power series is given by aout[][] = sum_i{al[i]*pl[i][][]} :

		// set initial values for the n entries in al[] and pal[] :
		fT wpf=1.0; //leading coefficient of power series
		for(i=0; i<n; ++i) {
			pal[i]=0;
			al[i]=wpf;
			wpf=opf(wpf,i+1); //compute (i+1)-th power series coefficent from the i-th coefficient, using the rule defined by opf()
			if(rescale==1) {
				//if matrix rescaling is used, next power series term will need additional factor of sfac compared to current term.
				wpf*=sfac;
			}
		}
		pal[n-1]=1.0;

		// next we iteratively add higher order power series terms to al[] till al[] stops changing
		// more precisely: the iteration will terminate after nhl_max consecutive iterations have not changed al[]	
		T ch,cho,tch; // temporary variables for iteratin
		int nhl=0; // counts the number of consecutive non-changing iterations  
		int nch; // counts the number of unchanged al[] coefficients in the curret iteration 
		for(j=n; j<mmax; ++j) {
			nch=0;
			ch=-pal[n-1]*crpl[0];
			cho=pal[0];
			pal[0]=ch;
			tch=al[0];
			al[0]+=wpf*ch;
			if(tch==al[0]) {
				++nch;
			}
			for(i=1; i<n; ++i) {
				ch=cho-pal[n-1]*crpl[i];
				cho=pal[i];
				pal[i]=ch;
				tch=al[i];
				al[i]+=wpf*ch;
				if(tch==al[i]) {
					++nch;
				}
			}
			if(nch>=n) {
				// no al[] coefficient has changed during current iteration
				++nhl;
				if(nhl>=nhl_max) {
					//terminate iteration
					break;
				}
			} else {
				// at least one al[] coefficent has changed during current iteration
				nhl=0;
			}

			wpf=opf(wpf,j+1); //compute (i+1)-th power series coefficent from i-th coefficient, using the rule defined by opf()

			if(rescale==1) {
				//if matrix rescaling is used, next power series term will need additional factor of sfac compared to current term.
				wpf*=sfac;
			}

		}

		// form output matrix by summing the 0-th to (n-1)-th matrix powers pl[] with corresponding weights al[] 
		for(i=0; i<n; ++i) {
			for(j=0; j<n; ++j) {
				aout[i][j]=al[0]*pl[0][i][j];
				for(k=1; k<n; ++k) {
					aout[i][j]+=al[k]*pl[k][i][j];
				}
			}
		}
	}

	void operator()(T** ain,T** dain,T** aout,T** daout,int rescale=0) {
		// applies the power series defined by opf() to the matrix ain[][], using the Cayley-Hamilton recursion, 
		// and writes the result to aout[][]; computes also the derivative of aout[][] with respect to ain[][] in
		// the direction of dain[][] and write the result to daout[][].

		int i,j,k;
		fT sfac=1.0; //scaling factor
		// compute the 0-th to n-th matrix powers of ta[][] :
		//   the i-th matrix power of ta[][] is stored in pl[i][][]
		//   the trace of the i-th matrix power of ta[][] is stored in trpl[i]
		set_to_identity(pl[0]);
		trpl[0]=n;
		if(rescale==1) {
			// if matrix rescaling is used, set sfac to be the magnitude of the largest element of ain[][]
			// and initiate the computation of the matrix powers from ain[][]/sfac instead of ain[][]:
			// (note that since we compute also the matrix powers pl[] from the rescaled matrix,
			// the Cayley-Hamilton coefficient a_{k,j}, with k=0,1,...,\infty; j=0,...,n-1 
			// will need to be rescaled only by a factor  sfac^{k}, instead of sfac^{k-j})
			sfac=max_matrix_el(ain);
			matrix_copy_scaled(ain,1.0/sfac,pl[1]);
		} else {
			matrix_copy(ain,pl[1]);
		}
		get_trace(pl[1],trpl[1]);
		for(i=2; i<=n; ++i) {
			j=i/2;
			k=i%2;
			matrix_mult_nn(pl[j],pl[j+k],pl[i],trpl[i]);
		}

		// compute the characteristic polynomial crpl[] from the traced powers trpl[]:
		crpl[n]=1;
		for(j=1; j<=n; ++j) {
			crpl[n-j]=0;
			for(i=1; i<=j; ++i) {
				crpl[n-j]-=crpl[n-(j-i)]*trpl[i];
			}
			crpl[n-j]/=j;
		}

		T*** tdpl=dpl[0];
		T* tdtrpl=dtrpl[0];
		T* tdcrpl=dcrpl[0];

		// compute the derivatives dpl[i][][] of the matrix powers pl[i][][] and corresponding traces dtrpl[i][][]
		set_to_zero(tdpl[0]);
		tdtrpl[0]=0;
		if(rescale==1) {
			//if we want to use wps factors from power siers also for the differentials,
			//therefore need to rescale the differential basis as well:
			matrix_copy_scaled(dain,1.0/sfac,tdpl[1]);
		} else {
			matrix_copy(dain,tdpl[1]);
		}
		get_trace(tdpl[1],tdtrpl[1]);
		if(n>=2) {
			matrix_mult_nn(pl[1],tdpl[1],tdpl[2]);  // dP_{2}=U.dP_{1}
			matrix_mult_nn(tdpl[1],pl[1],tmat); // tmat=dP_{1}.U
			matrix_add(tmat,tdpl[2],tdtrpl[2]); //dP_{2}+=tmat  -->  dP_{2}=U.dP_{1}+dP_{1}.U
			for(k=3; k<=n; ++k) {
				matrix_mult_nn(pl[1],tdpl[k-1],tdpl[k]); // dP_{k}=U.dP_{k-1}
				matrix_mult_nn_sub(pl[1],tmat,tdpl[k]); // dP_{k}-=U.tmat  -->  dP_{k}=U.dP_{k-1}-U.dP_{k-2}.U
				matrix_mult_nn(tdpl[k-1],pl[1],tmat);  // tmat=dP_{k-1}.U
				matrix_add(tmat,tdpl[k],tdtrpl[k]); //dP_{k}+=tmat  -->  dP_{k}=U.dP_{k-1}-U.dP_{k-2}.U+dP_{k-1}.U
			}
		}


		// compute the characteristic polynomial derivatives dcrpl[] 
		// from the traced powers trpl[] and their derivatives dtrpl[]:
		tdcrpl[n]=0;
		for(j=1; j<=n; ++j) {
			tdcrpl[n-j]=0;
			for(i=1; i<=j; ++i) {
				tdcrpl[n-j]-=crpl[n-(j-i)]*tdtrpl[i]+trpl[i]*tdcrpl[n-(j-i)];
			}
			tdcrpl[n-j]/=j;
		}

		// compute iteratively the n coefficients al[] and dal[] so that the Cayley-Hamilton result
		// for the matrix power series is given by aout[][] = sum_i{al[i]*pl[i][][]}
		// and its derivative by daout[][] = sum_i{dal[i]*pl[i][][]+al[i]*dpl[i][][]}:

		T* tdpal=dpal[0];
		T* tdal=dal[0];

		// set initial values for the n entries in al[] and pal[], as well as dal[] and dpal[] :
		fT wpf=1.0; //leading coefficient of power series 
		for(i=0; i<n; ++i) {
			pal[i]=0;
			al[i]=wpf;

			tdpal[i]=0;
			tdal[i]=0;
			wpf=opf(wpf,i+1); //compute (i+1)-th power series coefficent from the i-th coefficient, using the rule defined by opf()
			if(rescale==1) {
				//if matrix rescaling is used, next power series term will need additional factor of sfac compared to current term.
				wpf*=sfac;
			}
		}
		pal[n-1]=1.0;

		// next we iteratively add higher order power series terms to al[] and dal[] till al[] stops changing
		// more precisely: the iteration will terminate after nhl_max consecutive iterations have not changed al[]	
		T ch,cho,tch; // temporary variables for iteratin
		int nhl=0; // counts the number of consecutive non-changing iterations  
		int nch; // counts the number of unchanged al[] coefficients in the curret iteration 
		for(j=n; j<mmax; ++j) {
			// first update the derivative terms, since they depen on the current value of pal[n-1]:
			ch=-tdpal[n-1]*crpl[0]-pal[n-1]*tdcrpl[0];
			cho=tdpal[0];
			tdpal[0]=ch;
			tdal[0]+=wpf*ch;
			for(i=1; i<n; ++i) {
				ch=cho-tdpal[n-1]*crpl[i]-pal[n-1]*tdcrpl[i];
				cho=tdpal[i];
				tdpal[i]=ch;
				tdal[i]+=wpf*ch;
			}

			// now update the normal power series terms:
			nch=0;
			ch=-pal[n-1]*crpl[0];
			cho=pal[0];
			pal[0]=ch;
			tch=al[0];
			al[0]+=wpf*ch;
			if(tch==al[0]) {
				++nch;
			}
			for(i=1; i<n; ++i) {
				ch=cho-pal[n-1]*crpl[i];
				cho=pal[i];
				pal[i]=ch;
				tch=al[i];
				al[i]+=wpf*ch;
				if(tch==al[i]) {
					++nch;
				}
			}

			if(nch>=n) {
				// no al[] coefficient has changed during current iteration
				++nhl;
				if(nhl>=nhl_max) {
					//terminate iteration
					break;
				}
			} else {
				// at least one al[] coefficent has changed during current iteration
				nhl=0;
			}

			wpf=opf(wpf,j+1); //compute (i+1)-th power series coefficent from i-th coefficient, using the rule defined by opf()
			if(rescale==1) {
				//if matrix rescaling is used, next power series term will need additional factor of sfac compared to current term.
				wpf*=sfac;
			}
		}

		// form output matrices
		for(i=0; i<n; ++i) {
			for(j=0; j<n; ++j) {
				aout[i][j]=al[0]*pl[0][i][j];
				for(k=1; k<n; ++k) {
					aout[i][j]+=al[k]*pl[k][i][j];
				}
			}
		}

		for(i=0; i<n; ++i) {
			for(j=0; j<n; ++j) {
				daout[i][j]=tdal[0]*pl[0][i][j]+al[0]*tdpl[0][i][j];
				for(k=1; k<n; ++k) {
					daout[i][j]+=tdal[k]*pl[k][i][j]+al[k]*tdpl[k][i][j];
				}
			}
		}

	}

	void operator()(T** ain,T** aout,T**** daout,int rescale=0) {
		// applies the power series defined by opf() to the matrix ain[][], using the Cayley-Hamilton recursion, 
		// and writes the result to aout[][]; computes also the derivative of aout[][] with respect to each of 
		// the nxn components of ain[][] and write the result to daout[][][][] (the first two indices define the component
		// with respect to which the derivative is taken and the last two indices enumerate the components of the matrix-
		// valued derivative).

		int i,j,k;
		fT sfac=1.0; //scaling factor
		// compute the 0-th to n-th matrix powers of ta[][] :
		//   the i-th matrix power of ta[][] is stored in pl[i][][]
		//   the trace of the i-th matrix power of ta[][] is stored in trpl[i]
		set_to_identity(pl[0]);
		trpl[0]=n;
		if(rescale==1) {
			// if matrix rescaling is used, set sfac to be the magnitude of the largest element of ain[][]
			// and initiate the computation of the matrix powers from ain[][]/sfac instead of ain[][]:
			sfac=max_matrix_el(ain);
			matrix_copy_scaled(ain,1.0/sfac,pl[1]);
		} else {
			matrix_copy(ain,pl[1]);
		}
		get_trace(pl[1],trpl[1]);
		for(i=2; i<=n; ++i) {
			j=i/2;
			k=i%2;
			matrix_mult_nn(pl[j],pl[j+k],pl[i],trpl[i]);
		}

		// compute the characteristic polynomial crpl[] from the traced powers trpl[]:
		crpl[n]=1;
		for(j=1; j<=n; ++j) {
			crpl[n-j]=0;
			for(i=1; i<=j; ++i) {
				crpl[n-j]-=crpl[n-(j-i)]*trpl[i];
			}
			crpl[n-j]/=j;
		}

		T*** tdpl;
		T* tdtrpl;
		T* tdcrpl;
		int ic1,ic2;
		int idrv=0;
		for(ic1=0; ic1<n; ++ic1) {
			for(ic2=0; ic2<n; ++ic2) {
				// compute the derivatives dpl[i][][] of the matrix powers pl[i][][] and corresponding traces dtrpl[i][][]
				tdpl=dpl[idrv];
				tdtrpl=dtrpl[idrv];
				tdcrpl=dcrpl[idrv];

				set_to_zero(tdpl[0]);
				tdtrpl[0]=0;
				set_to_zero(tdpl[1]);
				if(rescale==1) {
					//if we want to use wps factors from power siers also for the differentials,
					//therefore need to rescale the differential basis as well:
					tdpl[1][ic1][ic2]=1.0/sfac;
					if(ic1==ic2) {
						tdtrpl[1]=1.0/sfac;
					} else {
						tdtrpl[1]=0;
					}
				} else {
					tdpl[1][ic1][ic2]=1.0;
					if(ic1==ic2) {
						tdtrpl[1]=1.0;
					} else {
						tdtrpl[1]=0;
					}
				}
				if(n>=2) {
					matrix_mult_nn(pl[1],tdpl[1],tdpl[2]);  // dP_{2}=U.dP_{1}
					matrix_mult_nn(tdpl[1],pl[1],tmat); // tmat=dP_{1}.U
					matrix_add(tmat,tdpl[2],tdtrpl[2]); //dP_{2}+=tmat  -->  dP_{2}=U.dP_{1}+dP_{1}.U
					for(k=3; k<=n; ++k) {
						matrix_mult_nn(pl[1],tdpl[k-1],tdpl[k]); // dP_{k}=U.dP_{k-1}
						matrix_mult_nn_sub(pl[1],tmat,tdpl[k]); // dP_{k}-=U.tmat  -->  dP_{k}=U.dP_{k-1}-U.dP_{k-2}.U
						matrix_mult_nn(tdpl[k-1],pl[1],tmat);  // tmat=dP_{k-1}.U
						matrix_add(tmat,tdpl[k],tdtrpl[k]); //dP_{k}+=tmat  -->  dP_{k}=U.dP_{k-1}-U.dP_{k-2}.U+dP_{k-1}.U
					}
				}


				// compute the characteristic polynomial derivatives dcrpl[] 
				// from the traced powers trpl[] and their derivatives dtrpl[]:
				tdcrpl[n]=0;
				for(j=1; j<=n; ++j) {
					tdcrpl[n-j]=0;
					for(i=1; i<=j; ++i) {
						tdcrpl[n-j]-=crpl[n-(j-i)]*tdtrpl[i]+trpl[i]*tdcrpl[n-(j-i)];
					}
					tdcrpl[n-j]/=j;
				}
				++idrv;
			}
		}

		// compute iteratively the n coefficients al[] and dal[] so that the Cayley-Hamilton result
		// for the matrix power series is given by aout[][] = sum_i{al[i]*pl[i][][]}
		// and its derivative by daout[][] = sum_i{dal[i]*pl[i][][]+al[i]*dpl[i][][]}:

		// set initial values for the n entries in al[] and pal[] :
		fT wpf=1.0; //leading coefficient of power series 
		for(i=0; i<n; ++i) {
			pal[i]=0;
			al[i]=wpf;
			wpf=opf(wpf,i+1); //compute (i+1)-th power series coefficent from the i-th coefficient, using the rule defined by opf()
			if(rescale==1) {
				//if matrix rescaling is used, next power series term will need additional factor of sfac compared to current term.
				wpf*=sfac;
			}
		}
		pal[n-1]=1.0;

		// set initial values for the entires in dal[][] and dpal[][] :
		T* tdpal;
		T* tdal;
		for(idrv=0; idrv<n*n; ++idrv) {
			tdpal=dpal[idrv];
			tdal=dal[idrv];
			for(i=0; i<n; ++i) {
				tdpal[i]=0;
				tdal[i]=0;
			}
		}


		// next we iteratively add higher order power series terms to al[] and the dal[][] till al[] stops changing
		// more precisely: the iteration will terminate after nhl_max consecutive iterations have not changed al[]	
		T ch,cho,tch; // temporary variables for iteration
		int nhl=0; // counts the number of consecutive non-changing iterations  
		int nch; // counts the number of unchanged al[] coefficients in the curret iteration 
		for(j=n; j<mmax; ++j) {

			// first update the derivative terms, since they depen on the current value of pal[n-1]:
			for(idrv=0; idrv<n*n; ++idrv) {
				tdcrpl=dcrpl[idrv];
				tdpal=dpal[idrv];
				tdal=dal[idrv];
				ch=-tdpal[n-1]*crpl[0]-pal[n-1]*tdcrpl[0];
				cho=tdpal[0];
				tdpal[0]=ch;
				tdal[0]+=wpf*ch;
				for(i=1; i<n; ++i) {
					ch=cho-tdpal[n-1]*crpl[i]-pal[n-1]*tdcrpl[i];
					cho=tdpal[i];
					tdpal[i]=ch;
					tdal[i]+=wpf*ch;
				}
			}

			// now update the normal power series terms:
			nch=0;
			ch=-pal[n-1]*crpl[0];
			cho=pal[0];
			pal[0]=ch;
			tch=al[0];
			al[0]+=wpf*ch;
			if(tch==al[0]) {
				++nch;
			}
			for(i=1; i<n; ++i) {
				ch=cho-pal[n-1]*crpl[i];
				cho=pal[i];
				pal[i]=ch;
				tch=al[i];
				al[i]+=wpf*ch;
				if(tch==al[i]) {
					++nch;
				}
			}

			if(nch>=n) {
				// no al[] coefficient has changed during current iteration
				++nhl;
				if(nhl>=nhl_max) {
					//terminate iteration
					break;
				}
			} else {
				// at least one al[] coefficent has changed during current iteration
				nhl=0;
			}

			wpf=opf(wpf,j+1); //compute (i+1)-th power series coefficent from i-th coefficient, using the rule defined by opf()
			if(rescale==1) {
				//if matrix rescaling is used, next power series term will need additional factor of sfac compared to current term.
				wpf*=sfac;
			}
		}

		// form output matrices
		for(i=0; i<n; ++i) {
			for(j=0; j<n; ++j) {
				aout[i][j]=al[0]*pl[0][i][j];
				for(k=1; k<n; ++k) {
					aout[i][j]+=al[k]*pl[k][i][j];
				}
			}
		}

		idrv=0;
		for(ic1=0; ic1<n; ++ic1) {
			for(ic2=0; ic2<n; ++ic2) {
				tdpl=dpl[idrv];
				tdal=dal[idrv];
				for(i=0; i<n; ++i) {
					for(j=0; j<n; ++j) {
						daout[ic1][ic2][i][j]=tdal[0]*pl[0][i][j]+al[0]*tdpl[0][i][j];
						for(k=1; k<n; ++k) {
							daout[ic1][ic2][i][j]+=tdal[k]*pl[k][i][j]+al[k]*tdpl[k][i][j];
						}
					}
				}
				++idrv;
			}
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

	fT max_matrix_el(T** lin1) {
		// find and return the component of lin1[][] of largest magnitutde
		int ic1,ic2;
		T* lin10;
		fT tres;
		fT res=0;
		for(ic1=0; ic1<n; ++ic1) {
			lin10=lin1[ic1];
			for(ic2=0; ic2<n; ++ic2) {
				tres=std::abs(lin10[ic2]);
				if(tres>res) {
					res=tres;
				}
			}
		}
		return res;
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
				lout0[ic2]=scalef*lin10[ic2];
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

private:
	int n;
	int mmax;
	int nhl_max;
	T** a;
	T*** pl;
	T**** dpl;
	T* trpl;
	T** dtrpl;
	T* crpl;
	T** dcrpl;
	T* pal;
	T** dpal;
	T* al;
	T** dal;
	T** tmat;
	fT* wps;
	fT(*opf)(fT,int);
};

