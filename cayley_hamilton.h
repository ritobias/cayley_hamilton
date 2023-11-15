#pragma once
#include <algorithm>
#include <complex>
#include <limits>
#include <cmath>

typedef double ftype;
typedef std::complex<ftype> ctype;

static const ftype _fmin=std::numeric_limits<ftype>::min();
static const ftype _fmax=std::numeric_limits<ftype>::max();
static const ftype _fprec=std::numeric_limits<ftype>::epsilon();
//static const ftype _lfprec=std::numeric_limits<ftype>::min_exponent;
static const ftype _lfprec=std::log(_fprec);

template<class T>
class cayley_hamilton {
public:
	cayley_hamilton(): n(0),a(0),pl(0),trpl(0),crpl(0),pal(0),tal(0),mmax(0) {
		opf=[](ftype pref,ftype x) { return pref/x; };
	}

	cayley_hamilton(int tn) {
		if(tn>0) {
			n=tn;
			new_matrix(a);
			new_matrix_array(pl,n+1);
			trpl=new T[n+1]();
			crpl=new T[n+1]();
			pal=new T[n]();
			tal=new T[n]();
			opf=[](ftype pref,ftype x) { return pref/x; };
			mmax=100*n;
		} else {
			n=0;
			a=0;
			pl=0;
			trpl=0;
			crpl=0;
			pal=new T[n]();
			tal=new T[n]();
			opf=[](ftype pref,ftype x) { return pref/x; };
			mmax=100;
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
		if(tal!=0) {
			delete[] tal;
			tal=0;
		}
		n=0;
	}

	void set_n(int tn) {
		if(tn>0) {
			if(tn!=n) {
				n=tn;
				if(a!=0) {
					delete_matrix(a);
				}
				new_matrix(a);
				if(pl!=0) {
					delete_matrix_array(pl,n+1);
				}
				new_matrix_array(pl,n+1);
				if(trpl!=0) {
					delete[] trpl;
					trpl=0;
				}
				trpl=new T[n+1];
				if(crpl!=0) {
					delete[] crpl;
					crpl=0;
				}
				crpl=new T[n+1];
				if(pal!=0) {
					delete[] pal;
					pal=0;
				}
				pal=new T[n];
				if(tal!=0) {
					delete[] tal;
					tal=0;
				}
				tal=new T[n];

				mmax=100*n;
			}
		} else {
			n=0;
			if(a!=0) {
				delete_matrix(a);
			}
			if(pl!=0) {
				delete_matrix_array(pl);
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
			if(tal!=0) {
				delete[] tal;
				tal=0;
			}
			mmax=0;
		}
	}

	void operator()(T** tain,T** taout) {
		int i,j,k;
		// compute the 0-th to n-th matrix powers of ta[][]
		// the i-th matrix powers of ta[][] is stored in pl[i][][]
		// the trace of the i-th matrix power of ta[][] is stored in trpl[i]
		set_to_identity(pl[0]);
		trpl[0]=n;
		matrix_copy(tain,pl[1]);
		get_trace(pl[1],trpl[1]);
		for(i=2; i<=n; ++i) {
			matrix_mult_nn(tain,pl[i-1],pl[i],trpl[i]);
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

		// compute iteratively the n coefficients tal[] so that the Cayley-Hamilton result
		// for the matrix power series is given by by: taout[][] = sum_i tal[i]*pl[i][][]

		//set initial values for the n entries in tal[]:
		ftype wpf=1.0; //leading coefficient of power series 
		for(i=0; i<n; ++i) {
			pal[i]=0;
			tal[i]=wpf;
			wpf=opf(wpf,(ftype)(i+1)); //compute (i+1)-th power series coefficent from i-th coefficient
		}
		pal[n-1]=1.0;

		T ch;
		T cho;
		T tch;
		int nhl_max=3;
		int nhl=nhl_max;
		int has_chn=0;

		// iteratively add higher order power series terms to tal[] till tal[] doesn't change anymore: 
		for(j=n; j<mmax; ++j) {
			has_chn=1;
			ch=-pal[n-1]*crpl[0];
			cho=pal[0];
			pal[0]=ch;
			tch=tal[0];
			tal[0]+=wpf*ch;
			if(tch==tal[0]) {
				has_chn=0;
			}
			for(i=1; i<n; ++i) {
				ch=cho-pal[n-1]*crpl[i];
				cho=pal[i];
				pal[i]=ch;
				tch=tal[i];
				tal[i]+=wpf*ch;
				if(tch==tal[i]) {
					has_chn=0;
				}
			}
			if(has_chn==0) {
				--nhl;
				if(nhl<=0) {
					break;
				}
			} else {
				nhl=nhl_max;
			}

			wpf=opf(wpf,(ftype)(j+1));
		}

		// form the output matrix:
		for(i=0; i<n; ++i) {
			for(j=0; j<n; ++j) {
				taout[i][j]=tal[0]*pl[0][i][j];
				for(k=1; k<n; ++k) {
					taout[i][j]+=tal[k]*pl[k][i][j];
				}
			}
		}
	}

private:
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

	void get_trace(T** lin1, T& rout) {
		rout=0;
		for(int ic1=0; ic1<n; ++ic1) {
			rout+=lin1[ic1][ic1];
		}
	}

	int n;
	int mmax;
	T** a;
	T*** pl;
	T* trpl;
	T* crpl;
	T* pal;
	T* tal;
	ftype(*opf)(ftype,ftype);
};

