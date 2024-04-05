#pragma once
#include<iostream>
#include <algorithm>
#include <complex>
#include <limits>
#include <cmath>

typedef double ftype;
typedef std::complex<ftype> ctype;

static const ftype _fprec=std::numeric_limits<ftype>::epsilon();

template<class T>
class cayley_hamilton {
public:
	cayley_hamilton(): n(0),a(0),pl(0),trpl(0),crpl(0),pal(0),al(0),mmax(0),nhl_max(0) {
		opf=[](ftype pref,int i) { return pref/(ftype)i; }; // returns the i-th coeff. of the power seris, computed from the (i-1)-th coeff. "pref" and "i"
													        // default rule generates the coeffs. for the powerseries of exp(), i.e. 1/i!
	}

	cayley_hamilton(int tn) {
		if(tn>0) {
			n=tn;
			new_matrix(a);
			new_matrix_array(pl,n+1);
			trpl=new T[n+1]();
			crpl=new T[n+1]();
			pal=new T[n]();
			al=new T[n]();
			opf=[](ftype pref,int i) { return pref/(ftype)i; };
			mmax=100*n;
			nhl_max=3;
		} else {
			n=0;
			a=0;
			pl=0;
			trpl=0;
			crpl=0;
			pal=0;
			al=0;
			opf=[](ftype pref,int i) { return pref/(ftype)i; };
			mmax=100;
			nhl_max=3;
		}
	}

	~cayley_hamilton() {
		if(a!=0) {
			delete_matrix(a);
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
		n=0;
	}

	void set_n(int tn) {
		//if an instance has been created with the trivial/empty constructor, 
		//this function has to be run afterwards to initialize the instance properly
		if(tn>0) {
			if(tn!=n) {
				if(a!=0) {
					delete_matrix(a);
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

				n=tn;

				new_matrix(a);
				if(pl!=0) {
					delete_matrix_array(pl,n+1);
				}
				new_matrix_array(pl,n+1);

				trpl=new T[n+1];

				crpl=new T[n+1];

				pal=new T[n];

				al=new T[n];

				mmax=100*n;
			}
		} else {
			if(a!=0) {
				delete_matrix(a);
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
			mmax=0;

			n=0;
		}
	}

	void operator()(T** ain,T** aout) {
		// applies the power series defined by opf() to the matrix ain[][], using the Cayley-Hamilton recursion, and writes the result to aout[][]
		int i,j,k;
		// compute the 0-th to n-th matrix powers of ta[][]
		// the i-th matrix powers of ta[][] is stored in pl[i][][]
		// the trace of the i-th matrix power of ta[][] is stored in trpl[i]
		set_to_identity(pl[0]);
		trpl[0]=n;
		matrix_copy(ain,pl[1]);
		get_trace(pl[1],trpl[1]);
		for(i=2; i<=n; ++i) {
			j=i/2;
			k=i%2;
			matrix_mult_nn(pl[j],pl[j+k],pl[i],trpl[i]);
		}

		// compute the characteristic polynomial crpl[] from the traced powers trpl[]
		crpl[n]=1;
		for(j=1; j<=n; ++j) {
			crpl[n-j]=0;
			for(i=1; i<=j; ++i) {
				crpl[n-j]-=crpl[n-(j-i)]*trpl[i];
			}
			crpl[n-j]/=j;
		}

		// compute iteratively the n coefficients al[] so that the Cayley-Hamilton result
		// for the matrix power series is given by: aout[][] = sum_i al[i]*pl[i][][]

		//set initial values for the n entries in al[] and pal[]:
		ftype wpf=1.0; //leading coefficient of power series 
		for(i=0; i<n; ++i) {
			pal[i]=0;
			al[i]=wpf;
			wpf=opf(wpf,i+1); //compute (i+1)-th power series coefficent from i-th coefficient, using the rule defined by opf()
		}
		pal[n-1]=1.0;

		// next we iteratively add higher order power series terms to al[] till al[] doesn't change anymore
		// more precisely: the iteration will terminate if nhl_max consecutive iterations have not changed al[]	
		T ch,cho,tch; // temporary storage for iteratin
		int nhl=0; // counts the number of consecutive non-changing iterations  
		int nch; // counts number of unchanged al[] coefficients in the curret iteration 
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
				// one or more al[] coefficent has changed during current iteration
				nhl=0;
			}

			wpf=opf(wpf,j+1); //compute (i+1)-th power series coefficent from i-th coefficient, using the rule defined by opf()
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

	void set_to_zero(T** ta) {
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

	void matrix_copy(T** lin1,T** lout) {
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

	void matrix_copy_a(T** lin1,T** lout) {
		int ic1,ic2;
		T* lin10;
		for(ic1=0; ic1<n; ++ic1) {
			lin10=lin1[ic1];
			for(ic2=0; ic2<n; ++ic2) {
				lout[ic2][ic1]=std::conj(lin10[ic2]);
			}
		}
	}

	void matrix_mult_nn(T** lin1,T** lin2,T** lout) {
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

	void get_trace(T** lin1,T& rout) {
		rout=0;
		for(int ic1=0; ic1<n; ++ic1) {
			rout+=lin1[ic1][ic1];
		}
	}

	void new_matrix(T**& ta,int init=-1) {
		ta=new T*[n];
		if(init==-1) {
			for(int ic1=0; ic1<n; ++ic1) {
				ta[ic1]=new T[n];
			}
		} else if(init==0) {
			for(int ic1=0; ic1<n; ++ic1) {
				ta[ic1]=new T[n]();
			}
		} else {
			for(int ic1=0; ic1<n; ++ic1) {
				ta[ic1]=new T[n]();
				ta[ic1][ic1]=1;
			}
		}
	}

	void delete_matrix(T**& ta) {
		for(int ic1=0; ic1<n; ++ic1) {
			delete[] ta[ic1];
		}
		delete[] ta;
		ta=0;
	}

	void new_matrix_array(T***& tnna,int len,int init=-1) {
		tnna=new T**[len];
		if(init==-1) {
			for(int i=0; i<len; ++i) {
				tnna[i]=new T*[n];
				for(int ic1=0; ic1<n; ++ic1) {
					tnna[i][ic1]=new T[n];
				}
			}
		} else if(init==0) {
			int ic1;
			for(int i=0; i<len; ++i) {
				tnna[i]=new T*[n];
				for(ic1=0; ic1<n; ++ic1) {
					tnna[i][ic1]=new T[n]();
				}
			}
		} else {
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
	T* trpl;
	T* crpl;
	T* pal;
	T* al;
	ftype(*opf)(ftype,int);
};


template<class T>
class sa_entry {
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
public:
	sun_algebra(int tn): ch(tn),tgen(0),telem0(0),telem(0) {
		n=tn;
		ngen=n*n-1;
		generator=new sparse_mat<ctype>[ngen];
		int i,j;
		for(i=0; i<ngen; ++i) {
			generator[i].init(n);
		}

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
				tgen.push(j,j,std::sqrt(2.0/((ftype)(i+1)*(i+2))));		
			}
			tgen.push(i+1,i+1,-std::sqrt(2.0*(ftype)(i+1)/(ftype)(i+2)));
			++igen;
		}
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

	void log_ah(ctype** inmat,ftype* outvec) {
		int i;
		for(i=0; i<ngen; ++i) {
			outvec[i]=0;
		}
		ctype nfc=ctype(0.0,0.5);
		ftype tiv;
		ctype ttiv;
		int maxit=5*n;
		ch.matrix_copy(inmat,tmat);
		int it;
		ftype outvsq=0;
		ftype tivsq;
		for(it=0; it<maxit; ++it) {
			tivsq=0;
			for(int igen=0; igen<ngen; ++igen) {
				tgen=generator+igen;
				telem0=tgen->elem;
				tiv=0;
				for(i=0; i<tgen->nelem; ++i) {
					telem=telem0+i;
					tiv+=std::imag(telem->val*tmat[telem->ind2][telem->ind1]);
				}
				outvec[igen]+=tiv;
				tivsq+=tiv*tiv;
			}
			outvsq+=tivsq;
			if(tivsq<10000.0*_fprec*outvsq) {
				break;
			}
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
			ch(tmat,tmat2); //matrix exponential
			ch.matrix_mult_nn(inmat,tmat2,tmat);
		}
		std::cout<<"iterations: "<<it<<std::endl;
	}


	int n;
	int ngen;
	sparse_mat<ctype>* generator;
	cayley_hamilton<ctype> ch;
	ctype** tmat;
	ctype** tmat2;
	ftype* tvec;
	sparse_mat<ctype>* tgen;
	sa_entry<ctype>* telem0;
	sa_entry<ctype>* telem;
};

