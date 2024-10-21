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

typedef double ftype; // float type to be used
typedef long double lftype; // long version of ftype

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

	cayley_hamilton(): n(0),pl(0),trpl(0),crpl(0),pal(0),al(0),kho(0),mmax(0),nhl_max(0),tmat1(0),tmat2(0) {
		//default constructor; will require a call to set_n() to set the size of the square matrices on which the class will operate

		opf=[](fT pref,int i) { return pref/(fT)i; }; // returns the i-th coeff. of the power seris, computed from the (i-1)-th coeff. "pref" and "i"
													        // default rule generates the coeffs. for the powerseries of exp(), i.e. 1/i!
	}

	cayley_hamilton(int tn,fT(*topf)(fT,int)=0) {
		// constructor with matrix size n=tn as argument and optional argument for function pointer or lambda that defines power series coeffients. 
		if(topf!=0) {
			// non-zero function pointer
			opf=topf;
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
			kho=new T[n]();
			opf=[](fT pref,int i) { return pref/(fT)i; };
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
			kho=0;
			opf=[](fT pref,int i) { return pref/(fT)i; };
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
		if(kho!=0) {
			delete[] kho;
			kho=0;
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
				if(kho!=0) {
					delete[] kho;
					kho=0;
				}

				n=tn;

				new_matrix(tmat1);

				new_matrix(tmat2);

				new_matrix_array(pl,n+1);

				trpl=new T[n+1];

				crpl=new T[n+1];

				pal=new T[n];

				al=new T[n];

				kho=new T[n];
				
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
			if(kho!=0) {
				delete[] kho;
				kho=0;
			}
			mmax=0;

			n=0;
		}
	}

	void set_opf(fT(*topf)(fT,int)) {
		//set opf to point to a user-defined function or lambda for generating the power series coefficients
		opf=topf;
	}

	/*
	void operator()(T** ain,T** aout,fT rescale=0) {
		// applies the power series defined by opf() to the matrix ain[][], using the Cayley-Hamilton recursion, 
		// and writes the result to aout[][]
		// if rescale is set to 1, then the computation will be performed with rescaled input matrix, which is
		// useful for matrices of large norm; if rescale isset to 0, no matrix rescaling will be performed. 
		if(n>0) {
			int i,j,k;
			fT sfac=1.0; //scaling factor
			// compute the 0-th to n-th matrix powers of ta[][] :
			//   the i-th matrix power of ta[][] is stored in pl[i][][]
			//   the trace of the i-th matrix power of ta[][] is stored in trpl[i]
			set_to_identity(pl[0]);
			trpl[0]=n;
			if(rescale>0) {
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
			for(i=2; i<=n; ++i) {
				j=i/2;
				k=i%2;
				matrix_mult_nn(pl[j],pl[j+k],pl[i],trpl[i]);
			}

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
			fT wpf=1.0; //leading coefficient of power series
			for(i=0; i<n; ++i) {
				pal[i]=0;
				al[i]=wpf;
				wpf=opf(wpf,i+1); //compute (i+1)-th power series coefficent from the i-th coefficient, using the rule defined by opf()
				if(rescale>0) {
					//if matrix rescaling is used, next power series term will need additional factor of sfac compared to current term.
					wpf*=sfac;
				}
			}
			pal[n-1]=1.0;

			// next we iteratively add higher order power series terms to al[] till al[] stops changing
			// more precisely: the iteration will terminate after nhl_max consecutive iterations have not changed al[]	
			T ch,cho,tch; // temporary variables for iterating
			int nhl=0; // counts the number of consecutive non-changing iterations  
			int nch; // counts the number of unchanged al[] coefficients in the curret iteration 
			fT s,rs=1.0; // used for normalizing the pal[]
			for(j=n; j<mmax; ++j) {
				nch=0;
				pal[n-1]*=rs;
				ch=-pal[n-1]*crpl[0];
				cho=pal[0]*rs;
				pal[0]=ch;
				s=std::norm(ch);
				tch=al[0];
				al[0]+=wpf*ch;
				if(tch==al[0]) {
					++nch;
				}
				for(i=1; i<n; ++i) {
					ch=cho-pal[n-1]*crpl[i];
					cho=pal[i]*rs;
					pal[i]=ch;
					s+=std::norm(ch);
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

				if(s>1.0) {
					// if s is bigger than 1, normalize pal[] by a factor rs=1.0/s in next itaration, and multiply wpf by s to compensate
					s=std::sqrt(s);
					wpf*=s;
					rs=1.0/s;
				} else {
					rs=1.0;
				}

				if(rescale>0) {
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
	}
	*/

	void operator()(T** ain,T** aout,fT rescale=0) {
		// applies the power series defined by opf() to the matrix ain[][], using the Cayley-Hamilton recursion, 
		// and writes the result to aout[][]
		// if rescale is set to 1, then the computation will be performed with rescaled input matrix, which is
		// useful for matrices of large norm; if rescale isset to 0, no matrix rescaling will be performed. 
		if(n>0) {
			int i,j,k;
			fT sfac=1.0; //scaling factor
			// compute the 0-th to n-th matrix powers of ta[][] :
			//   the i-th matrix power of ta[][] is stored in pl[i][][]
			//   the trace of the i-th matrix power of ta[][] is stored in trpl[i]
			set_to_identity(pl[0]);
			trpl[0]=n;
			if(rescale>0) {
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
			for(i=2; i<=n; ++i) {
				j=i/2;
				k=i%2;
				matrix_mult_nn(pl[j],pl[j+k],pl[i],trpl[i]);
			}

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
			fT wpf=1.0,twpf=1.0, ttwpf; //leading coefficient of power series and of running sum used for convergence check
			for(i=0; i<n; ++i) {
				pal[i]=0;
				al[i]=wpf;
				wpf=opf(wpf,i+1); //compute (i+1)-th power series coefficent from the i-th coefficient, using the rule defined by opf()
				if(rescale>0) {
					//if matrix rescaling is used, next power series term will need additional factor of sfac compared to current term.
					wpf*=sfac;
				}
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

				wpf=opf(wpf,j+1); //compute (i+1)-th power series coefficent from i-th coefficient, using the rule defined by opf()

				if(s>1.0) {
					// if s is bigger than 1, normalize pal[] by a factor rs=1.0/s in next itaration, and multiply wpf by s to compensate
					s=std::sqrt(s);
					wpf*=s;
					rs=1.0/s;
				} else {
					rs=1.0;
				}

				if(rescale>0) {
					//if matrix rescaling is used, next power series term will need additional factor of sfac compared to current term.
					wpf*=sfac;
				}

				ttwpf=twpf;
				twpf+=wpf;
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
	}


	void operator()(T** ain,T** aout,T** daout,fT rescale=0) {
		// applies the power series defined by opf() to the matrix ain[][], using the Cayley-Hamilton recursion, 
		// and writes the result to aout[][]; computes also the derivative of aout[][] in the direction of
		// daout[][] and overwrite daout[][] with the result :
		if(n>0) {
			// use unused entries of pl[] as storage for computation of differentaial components:
			T** kmats=pl[0];
			T** kh=pl[n];

			int i,j,k;
			fT sfac=1.0; //scaling factor
			// compute the 0-th to n-th matrix powers of ta[][] :
			//   the i-th matrix power of ta[][] is stored in pl[i][][]
			//   the trace of the i-th matrix power of ta[][] is stored in trpl[i]
			trpl[0]=n;
			if(rescale>0) {
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
			fT wpf=1.0,twpf=1.0,ttwpf,wpff; //leading coefficient of power series and of running sum used for convergence check
			set_to_zero(kmats);
			for(i=0; i<n; ++i) {
				pal[i]=0;
				al[i]=wpf;
				wpf=opf(wpf,i+1); //compute (i+1)-th power series coefficent from the i-th coefficient, using the rule defined by opf()

				k=i/2;
				for(j=i-k; j<=i; ++j) {
					kmats[i-j][j]=wpf;
				}

				if(rescale>0) {
					//if matrix rescaling is used, next power series term will need additional factor of sfac compared to current term.
					wpf*=sfac;
				}
			}
			pal[n-1]=1.0;
			// set initial values for the entries in kh[][] :
			set_to_zero(kh);
			k=(n-1)/2;
			wpff=opf(wpf,n+1);

			for(i=n-1-k; i<n; ++i) {
				kh[n-1-i][i]=wpff;
			}


			// next we iteratively add higher order power series terms to al[] till al[] stops changing
			// more precisely: the iteration will terminate after nhl_max consecutive iterations have not changed al[]	
			T cho; // temporary variables for iterating
			int nhl=0; // counts the number of consecutive non-changing iterations  
			fT s,rs=1.0; // used for normalizing the pal[]
			for(j=n; j<mmax; ++j) {
				s=0;
				cho=pal[n-1];
				for(i=n-1; i>0; --i) {
					pal[i]=pal[i-1]-cho*crpl[i];
					s+=std::norm(pal[i]);
					al[i]+=wpf*pal[i];
				}
				pal[0]=-cho*crpl[0];
				s+=std::norm(pal[0]);
				al[0]+=wpf*pal[0];

				wpf=opf(wpf,j+1); //compute (i+1)-th power series coefficent from i-th coefficient, using the rule defined by opf()

				if(s>1.0) {
					// if s is bigger than 1, normalize pal[] by a factor rs=1.0/s in next itaration, and multiply wpf by s to compensate
					s=std::sqrt(s);
					wpf*=s;
					rs=1.0/s;
				} else {
					rs=1.0;
				}

				for(i=0; i<n; ++i) {
					pal[i]*=rs;
				}

				wpff=opf(wpf,j+2)/wpf;
				if(rescale>0) {
					//if matrix rescaling is used, next power series term will need additional factor of sfac compared to current term.
					wpff*=sfac;
				}

				for(i=n-1; i>=0; --i) {
					cho=kh[i][n-1];
					for(k=n-1; k>i; --k) {
						kh[i][k]=kh[i][k-1]-cho*crpl[k];
						kmats[i][k]+=kh[i][k];
						kh[i][k]*=wpff;
					}
					if(i>0) {
						kh[i][i]=kh[i-1][i]-cho*crpl[i];
						kmats[i][i]+=kh[i][i];
						kh[i][i]*=wpff;
					} else {
						kh[0][0]=wpf*pal[0]-cho*crpl[0];
						kmats[0][0]+=kh[0][0];
						kh[0][0]*=wpff;
					}
				}

				if(rescale>0) {
					//if matrix rescaling is used, next power series term will need additional factor of sfac compared to current term.
					wpf*=sfac;
				}

				ttwpf=twpf;
				twpf+=wpf;
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
	}

	void operator()(T** ain,T** aout,T**** daout,fT rescale=0) {
		// applies the power series defined by opf() to the matrix ain[][], using the Cayley-Hamilton recursion, 
		// and writes the result to aout[][]; computes also the derivative of aout[][] with respect to each of 
		// the nxn components of ain[][] and write the result to daout[][][][] (the first two indices define the component
		// with respect to which the derivative is taken and the last two indices enumerate the components of the matrix-
		// valued derivative).
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
			if(rescale>0) {
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
			fT wpf=1.0, wpff; //leading coefficient of power series
			set_to_zero(kmats);
			for(i=0; i<n; ++i) {
				pal[i]=0;
				al[i]=wpf;
				wpf=opf(wpf,i+1); //compute (i+1)-th power series coefficent from the i-th coefficient, using the rule defined by opf()

				k=i/2;
				for(j=i-k; j<=i; ++j) {
					kmats[i-j][j]=wpf;
				}

				if(rescale>0) {
					//if matrix rescaling is used, next power series term will need additional factor of sfac compared to current term.
					wpf*=sfac;
				}
			}
			pal[n-1]=1.0;
			// set initial values for the entries in kh[][] :
			set_to_zero(kh);
			k=(n-1)/2;
			wpff=opf(wpf,n+1);
			for(i=n-1-k; i<n; ++i) {
				kh[n-1-i][i]=wpff;
			}


			// next we iteratively add higher order power series terms to al[] till al[] stops changing
			// more precisely: the iteration will terminate after nhl_max consecutive iterations have not changed al[]	
			T ch,cho,tch; // temporary variables for iterating
			int nhl=0; // counts the number of consecutive non-changing iterations  
			int nch; // counts the number of unchanged al[] coefficients in the curret iteration 
			fT s,rs=1.0; // used for normalizing the pal[]
			for(j=n; j<mmax; ++j) {
				nch=0;
				ch=-pal[n-1]*crpl[0];
				cho=pal[0];
				pal[0]=ch;
				s=std::norm(ch);
				tch=al[0];
				al[0]+=wpf*ch;
				if(tch==al[0]) {
					++nch;
				}
				for(i=1; i<n; ++i) {
					ch=cho-pal[n-1]*crpl[i];
					cho=pal[i];
					pal[i]=ch;
					s+=std::norm(ch);
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

				if(s>1.0) {
					// if s is bigger than 1, normalize pal[] by a factor rs=1.0/s in next itaration, and multiply wpf by s to compensate
					s=std::sqrt(s);
					wpf*=s;
					rs=1.0/s;
				} else {
					rs=1.0;
				}

				for(i=0; i<n; ++i) {
					pal[i]*=rs;
					kho[i]=wpf*pal[i];
				}

				wpff=opf(wpf,j+2)/wpf;
				if(rescale>0) {
					//if matrix rescaling is used, next power series term will need additional factor of sfac compared to current term.
					wpf*=sfac;
					wpff*=sfac;
				}
				for(i=0; i<n; ++i) {
					for(k=i; k<n; ++k) {
						ch=kho[k]-kh[k][n-1]*crpl[i];
						kho[k]=kh[i][k];
						kh[i][k]=ch*wpff;
						kmats[i][k]+=ch;
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
	}

	void operator()(T** ain,T* rout,T*** drout,fT rescale=0) {
		// applies the power series defined by opf() to the matrix ain[][], using the Cayley-Hamilton recursion, 
		// and writes Cayley-Hamilton coefficients to rout[]; computes also the derivatives of rout[] with respect
		// to the components of ain[][] and writes the result to the list of matrices drout[][][] :
		if(n>0) {
			// use unused entries of pl[] as storage for computation of differentaial components:
			T** kmats=pl[0];
			T** kh=pl[n];

			int i,j,k;
			fT sfac=1.0; //scaling factor
			// compute the 0-th to n-th matrix powers of ta[][] :
			//   the i-th matrix power of ta[][] is stored in pl[i][][]
			//   the trace of the i-th matrix power of ta[][] is stored in trpl[i]
			trpl[0]=n;
			if(rescale>0) {
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

			T*** dcrpl;
			new_matrix_array(dcrpl,n+1);
			set_to_identity_scaled(-crpl[n],dcrpl[n-1]);
			for(j=2; j<=n; ++j) {
				matrix_copy_scaled(pl[j-1],-(fT)j,dcrpl[n-j]);
				matrix_add(-crpl[n-(j-1)],dcrpl[n-j]);
				matrix_mult_scalar_add(dcrpl[n-(j-1)],trpl[1],dcrpl[n-j]);
				for(i=2; i<j; ++i) {
					matrix_mult_scalar_sub(pl[i-1],(fT)i*crpl[n-(j-i)],dcrpl[n-j]);
					matrix_mult_scalar_sub(dcrpl[n-(j-i)],trpl[i],dcrpl[n-j]);
				}
				matrix_mult_scalar(dcrpl[n-j],1.0/j,dcrpl[n-j]);
			}

			// compute iteratively the n coefficients al[] so that the Cayley-Hamilton result
			// for the matrix power series is given by aout[][] = sum_i{al[i]*pl[i][][]} :

			// set initial values for the n entries in al[] and pal[] and the nxn entries in kmats[][] :
			fT wpf=1.0; //leading coefficient of power series
			for(i=0; i<n; ++i) {
				pal[i]=0;
				rout[i]=wpf;
				wpf=opf(wpf,i+1); //compute (i+1)-th power series coefficent from the i-th coefficient, using the rule defined by opf()

				if(rescale>0) {
					//if matrix rescaling is used, next power series term will need additional factor of sfac compared to current term.
					wpf*=sfac;
				}
			}
			pal[n-1]=1.0;

			// set initial values for the entries in drout[][][] and dpal[][][]
			T*** dpal;
			new_matrix_array(dpal,n,0);

			// next we iteratively add higher order power series terms to al[] till al[] stops changing
			// more precisely: the iteration will terminate after nhl_max consecutive iterations have not changed al[]	
			T ch,cho,tch; // temporary variables for iterating
			T** dch;
			new_matrix(dch);
			T** dcho;
			new_matrix(dcho);
			int nhl=0; // counts the number of consecutive non-changing iterations  
			int nch; // counts the number of unchanged al[] coefficients in the curret iteration 
			fT s,rs=1.0; // used for normalizing the pal[]
			for(j=n; j<mmax; ++j) {
				pal[n-1]*=rs;
				matrix_mult_scalar(dpal[n-1],rs,dpal[n-1]);
				nch=0;
				ch=-pal[n-1]*crpl[0];

				matrix_mult_scalar(dpal[n-1],-crpl[0],dch);
				matrix_mult_scalar_sub(dcrpl[0],pal[n-1],dch);

				cho=pal[0]*rs;
				pal[0]=ch;

				matrix_mult_scalar(dpal[0],rs,dcho);
				matrix_copy(dch,dpal[0]);

				s=std::norm(ch);
				tch=rout[0];
				rout[0]+=wpf*ch;
				
				matrix_mult_scalar_add(dch,wpf,drout[0]);

				if(tch==rout[0]) {
					++nch;
				}
				for(i=1; i<n; ++i) {
					ch=cho-pal[n-1]*crpl[i];

					matrix_copy(dcho,dch);
					matrix_mult_scalar_sub(dpal[n-1],crpl[i],dch);
					matrix_mult_scalar_sub(dcrpl[i],pal[n-1],dch);

					cho=pal[i]*rs;
					pal[i]=ch;

					matrix_mult_scalar(dpal[i],rs,dcho);
					matrix_copy(dch,dpal[i]);

					s+=std::norm(ch);
					tch=rout[i];
					rout[i]+=wpf*ch;

					matrix_mult_scalar_add(dch,wpf,drout[i]);

					if(tch==rout[i]) {
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

				if(s>1.0) {
					// if s is bigger than 1, normalize pal[] by a factor rs=1.0/s in next itaration, and multiply wpf by s to compensate
					s=std::sqrt(s);
					wpf*=s;
					rs=1.0/s;
				} else {
					rs=1.0;
				}

				if(rescale>0) {
					//if matrix rescaling is used, next power series term will need additional factor of sfac compared to current term.
					wpf*=sfac;
				}

			}

			if(rescale>0) {
				wpf=rescale;
				matrix_mult_scalar(drout[0],wpf,drout[0]);
				for(k=1; k<n; ++k) {
					rout[k]*=wpf;
					wpf*=rescale;
					matrix_mult_scalar(drout[k],wpf,drout[k]);
				}
			}

			delete_matrix(dch);
			delete_matrix(dcho);
			delete_matrix_array(dpal,n);
			delete_matrix_array(dcrpl,n+1);
		}
	}

	void operator()(T** ain,T** aout,T*** drout,fT rescale=0) {
		// applies the power series defined by opf() to the matrix ain[][], using the Cayley-Hamilton recursion, 
		// and writes the result to aout[][]; computes also the derivative of aout[][] in the direction of
		// daout[][] and overwrite daout[][] with the result :
		if(n>0) {
			// use unused entries of pl[] as storage for computation of differentaial components:
			T** kmats=pl[0];
			T** kh=pl[n];

			int i,j,k;
			fT sfac=1.0; //scaling factor
			// compute the 0-th to n-th matrix powers of ta[][] :
			//   the i-th matrix power of ta[][] is stored in pl[i][][]
			//   the trace of the i-th matrix power of ta[][] is stored in trpl[i]
			trpl[0]=n;
			if(rescale>0) {
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
			fT wpf=1.0,wpff; //leading coefficient of power series
			set_to_zero(kmats);
			for(i=0; i<n; ++i) {
				pal[i]=0;
				al[i]=wpf;
				wpf=opf(wpf,i+1); //compute (i+1)-th power series coefficent from the i-th coefficient, using the rule defined by opf()

				k=i/2;
				for(j=i-k; j<=i; ++j) {
					kmats[i-j][j]=wpf;
				}

				if(rescale>0) {
					//if matrix rescaling is used, next power series term will need additional factor of sfac compared to current term.
					wpf*=sfac;
				}
			}
			pal[n-1]=1.0;
			// set initial values for the entries in kh[][] :
			set_to_zero(kh);
			k=(n-1)/2;
			wpff=opf(wpf,n+1);

			for(i=n-1-k; i<n; ++i) {
				kh[n-1-i][i]=wpff;
			}


			// next we iteratively add higher order power series terms to al[] till al[] stops changing
			// more precisely: the iteration will terminate after nhl_max consecutive iterations have not changed al[]	
			T ch,cho,tch; // temporary variables for iterating
			int nhl=0; // counts the number of consecutive non-changing iterations  
			int nch; // counts the number of unchanged al[] coefficients in the curret iteration 
			fT s,rs=1.0; // used for normalizing the pal[]
			for(j=n; j<mmax; ++j) {
				nch=0;
				ch=-pal[n-1]*crpl[0];
				cho=pal[0];
				pal[0]=ch;
				s=std::norm(ch);
				tch=al[0];
				al[0]+=wpf*ch;
				if(tch==al[0]) {
					++nch;
				}
				for(i=1; i<n; ++i) {
					ch=cho-pal[n-1]*crpl[i];
					cho=pal[i];
					pal[i]=ch;
					s+=std::norm(ch);
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

				if(s>1.0) {
					// if s is bigger than 1, normalize pal[] by a factor rs=1.0/s in next itaration, and multiply wpf by s to compensate
					s=std::sqrt(s);
					wpf*=s;
					rs=1.0/s;
				} else {
					rs=1.0;
				}

				for(i=0; i<n; ++i) {
					pal[i]*=rs;
					kho[i]=wpf*pal[i];
				}

				wpff=opf(wpf,j+2)/wpf;
				if(rescale>0) {
					//if matrix rescaling is used, next power series term will need additional factor of sfac compared to current term.
					wpf*=sfac;
					wpff*=sfac;
				}
				for(i=0; i<n; ++i) {
					for(k=i; k<n; ++k) {
						ch=kho[k]-kh[k][n-1]*crpl[i];
						kho[k]=kh[i][k];
						kh[i][k]=ch*wpff;
						kmats[i][k]+=ch;
					}
				}

			}


			// form output matrix by summing the 0-th to (n-1)-th matrix powers pl[] with corresponding weights al[] 
			set_to_identity_scaled(al[0],aout);
			for(j=1; j<n; ++j) {
				matrix_add_scaled(pl[j],al[j],aout);
			}

			set_to_zero(tmat1);
			for(i=0; i<n; ++i) {
				tmat1[i][0]=-kmats[i][n-1];
			}
			for(j=1; j<n; ++j) {
				cho=tmat1[n-1][j-1];
				for(i=n-1; i>0; --i) {
					tmat1[i][j]=tmat1[i-1][j-1]-cho*crpl[i];
				}
				tmat1[0][j]=-cho*crpl[0];
			}

			T*** dcrpl;
			new_matrix_array(dcrpl,n+1);
			set_to_identity_scaled(-crpl[n],dcrpl[n-1]);
			for(j=2; j<=n; ++j) {
				matrix_copy_scaled(pl[j-1],-(fT)j,dcrpl[n-j]);
				matrix_add(-crpl[n-(j-1)],dcrpl[n-j]);
				matrix_mult_scalar_add(dcrpl[n-(j-1)],trpl[1],dcrpl[n-j]);
				for(i=2; i<j; ++i) {
					matrix_mult_scalar_sub(pl[i-1],(fT)i*crpl[n-(j-i)],dcrpl[n-j]);
					matrix_mult_scalar_sub(dcrpl[n-(j-i)],trpl[i],dcrpl[n-j]);
				}
				matrix_mult_scalar(dcrpl[n-j],1.0/j,dcrpl[n-j]);
			}
			
			for(k=0; k<n; ++k) {
				matrix_mult_scalar(dcrpl[0],tmat1[k][0],drout[k]);
				for(j=1; j<n; ++j) {
					matrix_mult_scalar_add(dcrpl[j],tmat1[k][j],drout[k]);
				}
			}

			delete_matrix_array(dcrpl,n+1);



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
				lout0[ic2]=lin10[ic2]*scalef;
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
	T*** pl;
	T* trpl;
	T* crpl;
	T* pal;
	T* kho;
	T* al;
	T** tmat1;
	T** tmat2;
	fT(*opf)(fT,int);
};

