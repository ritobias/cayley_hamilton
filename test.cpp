/*
    cpp file for testing the sun_algebra and cayley_hamilton classes
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
#include<iostream>
#include<random>
#include <fstream>
#include <sstream>
#include "cayley_hamilton.h"
#include "sun_algebra.h"
#include <dirent.h>
#include <map>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/timeb.h>
#include <ctime>
#include <deque>
#include <algorithm>

int getmcount() {
    timeb tb;
    ftime(&tb);
    int nc=tb.millitm+(tb.time&0xfffff)*1000;
    return nc;
}

int getmspan(int nts) {
    int ns=getmcount()-nts;
    if(ns<0)
        ns+=0x100000*1000;
    return ns;
}

int main()
{
    // Run some tests on the sun_algebra and/or the cayley_hamilton class:
    std::random_device r;
    std::seed_seq seedsq{r(),r(),r(),r(),r(),r(),r(),r()};
    std::mt19937_64 rng(seedsq);
    rng.discard(10000-1);
    std::uniform_real_distribution<> dis(0,1.0);
    std::normal_distribution<> ndis{1.0,1.0};

    std::cout<<"_fprec="<<_fprec<<std::endl;

    if(0) {
        //run tests on the sun_algebra class
        int n=7; // matrix size

        //initialize an instance of the sun_algebra class with matrix size n (i.e. Lie algebra su(n))
        sun_algebra sun(n);

        //allocate some matrices using utility functions from the cayley_hamilton sub-class of sun_algebra:
        ctype** a;
        sun.ch.new_matrix(a);
        ctype** ea;
        sun.ch.new_matrix(ea);

        // specify a random vector "v" of magnitude 2*alpha in space spanned by Lie algebra basis
        ftype alpha=0.5; //value for which "v" points to point within injectivty domain of exponential map exp_{id} : su(n)=T_{id}SU(n) -> SU(n))
        //ftype alpha=5.0; //value for which "v" points to point outside of injectivity domain of exponetial map exp_{id} : su(n)=T_{id}SU(n) -> SU(n))
        ftype* v=new ftype[sun.ngen];
        ftype vnorm=0;
        for(int i=0; i<sun.ngen; ++i) {
            v[i]=(ftype)ndis(rng);
            vnorm+=v[i]*v[i];
        }
        vnorm=std::sqrt(vnorm);
        for(int i=0; i<sun.ngen; ++i) {
            v[i]*=(ftype)2.0*alpha/vnorm;
        }

        // get the anit-hermitian matrix representation of "v": a=i*v.T ,  where T={lambda_0 / 2, ..., lambda_ngen / 2}
        sun.get_alg_mat_ah(v,a);

        std::cout<<"a:"<<std::endl;
        sun.ch.print_matrix(a,7);
        std::cout<<std::endl;

        // get the group matrix ea=Exp(i*v.T) 
        sun.get_grp_mat(v,ea);

        std::cout<<"ea=exp(a):"<<std::endl;
        sun.ch.print_matrix(ea,7);
        std::cout<<std::endl;

        // take the log of ea by determining the vector v in Lie algebra space for which ea=Exp(i*v.T) 
        // (note that the result of this does not need to agree with the initial "v" if "alpha" has been chose sufficently large,
        // so that "v" points outside of the injectivty domain of the exponential map exp_{id} : su(n)=T_{id}SU(n) -> SU(n))
        sun.get_log(ea,v);

        // get the corresponding anti hermitian matrix a=i*v.T
        sun.get_alg_mat_ah(v,a);

        std::cout<<"lea=log(ea):"<<std::endl;
        sun.ch.print_matrix(a,7);
        std::cout<<std::endl;

        // get the group matrix ea=Exp(i*v.T)
        // (this should now agree again with the group matrix that was obtained with the original "v") 
        sun.get_grp_mat(v,ea);

        std::cout<<"ea=exp(lea):"<<std::endl;
        sun.ch.print_matrix(ea,7);
        std::cout<<std::endl;

        //clean up memory used for sun_algebra tests:
        sun.ch.delete_matrix(a);
        sun.ch.delete_matrix(ea);
        delete[] v;
    }

    if(0) {
        //run tests on the cayley_hamilton class
        int n=5; // matrix size

        //initialize an instance of the sun_algebra class with matrix size n (i.e. Lie algebra su(n))
        cayley_hamilton<ctype> ch(n);
        chexp<ctype> mexp(n);

        //allocate some matrices using utility functions from the cayley_hamilton sub-class of sun_algebra:
        ctype** a;
        ch.new_matrix(a);
        ctype** ea;
        ch.new_matrix(ea);

        //test cayley_hamilton with a specific matrix to compare with known results  
        a[0][0]=ctype(0,0.850177714654586580);
        a[0][1]=ctype(0.274160860652971849,0.032390014426889491);
        a[0][2]=ctype(0.340537026101288887,0.180479527980594190);
        a[0][3]=ctype(1.000000000000000000,0.369532334946793384);
        a[0][4]=ctype(0.1010476580952816006,0.0120220519818148297);
        a[1][0]=ctype(-0.274160860652971849,0.032390014426889491);
        a[1][1]=ctype(0,0.1073405743445747240);
        a[1][2]=ctype(0.240482580987094519,0.169218533199260390);
        a[1][3]=ctype(0.0793172650637539167,0.0564989168424858690);
        a[1][4]=ctype(0.484908776368710104,0.071918635695438216);
        a[2][0]=ctype(-0.340537026101288887,0.180479527980594190);
        a[2][1]=ctype(-0.240482580987094519,0.169218533199260390);
        a[2][2]=ctype(0,0.349934243921792953);
        a[2][3]=ctype(0.294899468664304942,0.292568842786627870);
        a[2][4]=ctype(0.164813462851289394,0.252827614108137921);
        a[3][0]=ctype(-1.000000000000000000,0.369532334946793384);
        a[3][1]=ctype(-0.0793172650637539167,0.0564989168424858690);
        a[3][2]=ctype(-0.294899468664304942,0.292568842786627870);
        a[3][3]=ctype(0,-0.788927360417937651);
        a[3][4]=ctype(0.493527673285568902,0.484649871927598008);
        a[4][0]=ctype(-0.1010476580952816006,0.0120220519818148297);
        a[4][1]=ctype(-0.484908776368710104,0.071918635695438216);
        a[4][2]=ctype(-0.164813462851289394,0.252827614108137921);
        a[4][3]=ctype(-0.493527673285568902,0.484649871927598008);
        a[4][4]=ctype(0,-0.518525172503016440);

        //rescale matrix by factor rsf:
        ftype rsf=1.0;
        ch.matrix_mult_scalar(a,rsf,a);

        std::cout<<"a:"<<std::endl;
        ch.print_matrix(a,7);

        // determine the matrix exponential ea[][] of a[][]:
        ch(a,ea);
        std::cout<<"ea=exp(a):"<<std::endl;
        ch.print_matrix(ea,7);

        // determine the derivative of the exponential ea[][] of a[][] in the direction of a[][] itself:
        ctype** tdea;
        ch.new_matrix(tdea);
        ch.matrix_copy(a,tdea);
        ch(a,ea,tdea,2.0/rsf);
        std::cout<<"ea=exp(a):"<<std::endl;
        ch.print_matrix(ea,7);

        std::cout<<"tdea=dexp(a)(a):"<<std::endl;
        ch.print_matrix(tdea,7);

        ctype** tdeab;
        ch.new_matrix(tdeab);
        ch.matrix_copy(a,tdeab);
        mexp(a,ea,tdeab);
        std::cout<<"eab=exp(a):"<<std::endl;
        ch.print_matrix(ea,7);

        std::cout<<"tdeab=dexp(a)(a):"<<std::endl;
        ch.print_matrix(tdeab,7);

        ch.delete_matrix(tdea);
        ch.delete_matrix(tdeab);


        if(1) {
            // determine the Cayley-Hamilton coefficients of the exp(a) and their differentials:
            ctype* ra=new ctype[n];
            ctype** ka;
            ch.new_matrix(ka);
            ch.get_r_k(a,ra,ka);

            ctype* ra2=new ctype[n];
            ctype** ka2;
            ch.new_matrix(ka2);
            mexp.get_r_k(a,ra2,ka2);

            std::cout<<FIXED_FLOAT(8);
            std::cout<<"ra:"<<std::endl;
            for(int i=0; i<n; ++i) {
                std::cout<<ra[i]<<" ";
            }
            std::cout<<std::endl;
            std::cout<<"ra2:"<<std::endl;
            for(int i=0; i<n; ++i) {
                std::cout<<ra2[i]<<" ";
            }
            std::cout<<std::endl;
            std::cout<<"ra-ra2:"<<std::endl;
            for(int i=0; i<n; ++i) {
                std::cout<<ra[i]-ra2[i]<<" ";
            }
            std::cout<<std::endl;
            std::cout<<"ka[][]:"<<std::endl;
            ch.print_matrix(ka,8);
            std::cout<<"ka2[][]:"<<std::endl;
            ch.print_matrix(ka2,8);
            for(int i=0; i<n; ++i) {
                for(int j=0; j<n; ++j) {
                    ka[i][j]-=ka2[i][j];
                }
            }
            std::cout<<"ka[][]-ka2[][]:"<<std::endl;
            ch.print_matrix(ka,8);

            delete[] ra;
            ch.delete_matrix(ka);
            delete[] ra2;
            ch.delete_matrix(ka2);
        }



        if(1) {
            // determine the Cayley-Hamilton coefficients of the exp(a) and their differentials:
            ctype* ra=new ctype[n];
            ctype*** dra;
            ch.new_matrix_array(dra,n);
            ch.get_r_dr(a,ra,dra);

            ctype* ra2=new ctype[n];
            ctype*** dra2;
            ch.new_matrix_array(dra2,n);
            mexp.get_r_dr(a,ra2,dra2);
            std::cout<<"ra:"<<std::endl;
            for(int i=0; i<n; ++i) {
                std::cout<<ra[i]<<" ";
            }
            std::cout<<std::endl;
            std::cout<<"ra2:"<<std::endl;
            for(int i=0; i<n; ++i) {
                std::cout<<ra2[i]<<" ";
            }
            std::cout<<std::endl;
            std::cout<<"ra-ra2:"<<std::endl;
            for(int i=0; i<n; ++i) {
                std::cout<<ra[i]-ra2[i]<<" ";
            }
            std::cout<<std::endl;

            for(int k=0; k<n; ++k) {
                std::cout<<"dra["<<k<<"]:"<<std::endl;
                ch.print_matrix(dra[k],8);
                std::cout<<"dra2["<<k<<"]:"<<std::endl;
                ch.print_matrix(dra2[k],8);
                std::cout<<"dra["<<k<<"]-dra2["<<k<<"]:"<<std::endl;
                ch.matrix_sub(dra2[k],dra[k]);
                ch.print_matrix(dra[k],8);
                std::cout<<std::endl;
            }

            if(0) {
                std::cout<<std::endl;
                std::cout<<"dr:"<<std::endl;
                ctype** tmat;
                ch.new_matrix(tmat);
                ctype** tapow;
                ch.new_matrix(tapow);
                ctype** ttapow;
                ch.new_matrix(ttapow);
                for(int ic1=0; ic1<n; ++ic1) {
                    for(int ic2=0; ic2<n; ++ic2) {
                        std::cout<<"dr["<<ic1<<"]["<<ic2<<"]:"<<std::endl;
                        ch.set_to_zero(tmat);
                        ch.set_to_identity(tapow);
                        for(int k=0; k<n; ++k) {
                            ch.matrix_mult_scalar_add(tapow,dra[k][ic1][ic2],tmat);
                            ch.matrix_mult_nn(tapow,a,ttapow);
                            ch.matrix_copy(ttapow,tapow);
                        }
                        ch.print_matrix(tmat,7);
                        std::cout<<std::endl;
                        std::cout<<std::endl;

                    }
                    std::cout<<std::endl;
                }
                std::cout<<std::endl;
                ch.delete_matrix(tmat);
                ch.delete_matrix(tapow);
                ch.delete_matrix(ttapow);
            }
            delete[] ra;
            ch.delete_matrix_array(dra,n);
            delete[] ra2;
            ch.delete_matrix_array(dra2,n);
        }

        if(1) {
            // allocate memory for the differential matrices
            ctype**** dea=new ctype***[n];
            for(int ic1=0; ic1<n; ++ic1) {
                dea[ic1]=new ctype**[n];
                for(int ic2=0; ic2<n; ++ic2) {
                    ch.new_matrix(dea[ic1][ic2]);
                }
            }

            // determine the matrix exponential ea[][] of a[][] and the corresponding set differentials dea[i][j][][]=dea[][]/da[i][j]:
            ch(a,ea,dea,0.5);

            // allocate more memory to compare the results obtained with and wouthout scaling of the input matrix (should yield the same results):
            ctype** ea2;
            ch.new_matrix(ea2);
            ctype**** dea2=new ctype***[n];
            for(int ic1=0; ic1<n; ++ic1) {
                dea2[ic1]=new ctype**[n];
                for(int ic2=0; ic2<n; ++ic2) {
                    ch.new_matrix(dea2[ic1][ic2]);
                }
            }

            // determine the matrix exponential ea2[][] of a[][] and the corresponding of set differentials dea2[i][j][][]=dea2[][]/da[j][i]
            // using input matrix rescaling:
            mexp(a,ea2,dea2);

            // print the results for the matrix exponentials:
            std::cout<<"ea=exp(a):"<<std::endl;
            ch.print_matrix(ea,7);

            std::cout<<"ea2=exp(a):"<<std::endl;
            ch.print_matrix(ea2,7);

            std::cout<<"ea-ea2:"<<std::endl;
            for(int i=0; i<n; ++i) {
                for(int j=0; j<n; ++j) {
                    ea[i][j]-=ea2[i][j];
                }
            }
            ch.print_matrix(ea,7);

            // print the results for the derivative terms:
            for(int ic1=0; ic1<n; ++ic1) {
                for(int ic2=0; ic2<n; ++ic2) {
                    std::cout<<"dea["<<ic1<<"]["<<ic2<<"]:"<<std::endl;
                    ch.print_matrix(dea[ic1][ic2],7);

                    std::cout<<"dea2["<<ic1<<"]["<<ic2<<"]:"<<std::endl;
                    ch.print_matrix(dea2[ic1][ic2],7);

                    std::cout<<"dea["<<ic1<<"]["<<ic2<<"]-dea2["<<ic1<<"]["<<ic2<<"]:"<<std::endl;
                    for(int i=0; i<n; ++i) {
                        for(int j=0; j<n; ++j) {
                            dea[ic1][ic2][i][j]-=dea2[ic1][ic2][i][j];
                        }
                    }
                    ch.print_matrix(dea[ic1][ic2],7);
                }
            }

            for(int ic1=0; ic1<n; ++ic1) {
                for(int ic2=0; ic2<n; ++ic2) {
                    ch.delete_matrix(dea[ic1][ic2]);
                }
                delete[] dea[ic1];
            }
            delete[] dea;

            ch.delete_matrix(ea2);

            for(int ic1=0; ic1<n; ++ic1) {
                for(int ic2=0; ic2<n; ++ic2) {
                    ch.delete_matrix(dea2[ic1][ic2]);
                }
                delete[] dea2[ic1];
            }
            delete[] dea2;
        }

        //clean up memory used for cayley_hamilton tests:
        ch.delete_matrix(a);
        ch.delete_matrix(ea);
    }

    if(0) {
        //run tests on the cayley_hamilton class
        int n=2; // matrix size

        //initialize an instance of the sun_algebra class with matrix size n (i.e. Lie algebra su(n))
        cayley_hamilton<ftype> ch(n);

        //allocate some matrices using utility functions from the cayley_hamilton sub-class of sun_algebra:
        ftype** a;
        ch.new_matrix(a);
        ftype** ea;
        ch.new_matrix(ea);

        a[0][0]=-49.0;
        a[0][1]=24.0;
        a[1][0]=-64.0;
        a[1][1]=31.0;


        // print out the matrix a[][]:
        std::cout<<"a:"<<std::endl;
        ch.print_matrix(a,7);

        // determine the matrix exponential ea[][] of a[][]:
        int niter=ch(a,ea);
        std::cout<<"ea=exp(a):"<<std::endl;
        ch.print_matrix(ea,7);
        std::cout<<"niter: "<<niter<<std::endl;
        std::cout<<std::endl;

        //clean up memory used for cayley_hamilton tests:
        ch.delete_matrix(a);
        ch.delete_matrix(ea);
    }


    if(0) {
        //run tests on the cayley_hamilton class
        int n=4; // matrix size

        //initialize an instance of the sun_algebra class with matrix size n (i.e. Lie algebra su(n))
        cayley_hamilton<ftype> ch(n);

        //allocate some matrices using utility functions from the cayley_hamilton sub-class of sun_algebra:
        ftype** a;
        ch.new_matrix(a,0);
        ftype** ea;
        ch.new_matrix(ea);

        a[0][1]=6.0;
        a[1][2]=6.0;
        a[2][3]=6.0;


        // print out the matrix a[][]:
        std::cout<<"a:"<<std::endl;
        ch.print_matrix(a,7);

        // determine the matrix exponential ea[][] of a[][]:
        int niter=ch(a,ea);
        std::cout<<"ea=exp(a):"<<std::endl;
        ch.print_matrix(ea,7);
        std::cout<<"niter: "<<niter<<std::endl;
        std::cout<<std::endl;

        //clean up memory used for cayley_hamilton tests:
        ch.delete_matrix(a);
        ch.delete_matrix(ea);
    }


    if(0) {
        //run tests on the cayley_hamilton class
        int n=2; // matrix size

        //initialize an instance of the sun_algebra class with matrix size n (i.e. Lie algebra su(n))
        cayley_hamilton<ftype> ch(n);

        //allocate some matrices using utility functions from the cayley_hamilton sub-class of sun_algebra:
        ftype** a;
        ch.new_matrix(a,0);
        ftype** ea;
        ch.new_matrix(ea);

        a[0][0]=1.0+1.0e-5;
        a[0][1]=1.0;
        a[1][1]=1.0-1.0e-5;

        // print out the matrix a[][]:
        std::cout<<"a:"<<std::endl;
        ch.print_matrix(a,7);

        // determine the matrix exponential ea[][] of a[][]:
        int niter=ch(a,ea);
        std::cout<<"ea=exp(a):"<<std::endl;
        ch.print_matrix(ea,7);
        std::cout<<"niter: "<<niter<<std::endl;
        std::cout<<std::endl;

        //clean up memory used for cayley_hamilton tests:
        ch.delete_matrix(a);
        ch.delete_matrix(ea);
    }

    //run tests on the cayley_hamilton class
    if(0) {
        std::cout<<"test #1:"<<std::endl;
        int n=2; // matrix size

        //initialize an instance of the sun_algebra class with matrix size n (i.e. Lie algebra su(n))
        chexp<ftype> ch(n);

        //allocate some matrices using utility functions from the cayley_hamilton sub-class of sun_algebra:
        ftype** a;
        ch.new_matrix(a,0);
        ftype** ea;
        ch.new_matrix(ea);

        a[0][0]=1.0;
        a[1][1]=2.0;

        // print out the matrix a[][]:
        std::cout<<"a:"<<std::endl;
        ch.print_matrix(a,7);

        // determine the matrix exponential ea[][] of a[][]:
        int niter=ch(a,ea);
        std::cout<<"ea=exp(a):"<<std::endl;
        ch.print_matrix(ea,7);
        std::cout<<"niter: "<<niter<<std::endl;
        std::cout<<std::endl;

        //clean up memory used for cayley_hamilton tests:
        ch.delete_matrix(a);
        ch.delete_matrix(ea);
    }

    if(0) {
        std::cout<<"test #2:"<<std::endl;
        int n=2; // matrix size

        //initialize an instance of the sun_algebra class with matrix size n (i.e. Lie algebra su(n))
        chexp<ftype> ch(n);

        //allocate some matrices using utility functions from the cayley_hamilton sub-class of sun_algebra:
        ftype** a;
        ch.new_matrix(a,0);
        ftype** ea;
        ch.new_matrix(ea);

        a[0][0]=1.0;
        a[0][1]=3.0;
        a[1][0]=3.0;
        a[1][1]=2.0;


        // print out the matrix a[][]:
        std::cout<<"a:"<<std::endl;
        ch.print_matrix(a,7);

        // determine the matrix exponential ea[][] of a[][]:
        int niter=ch(a,ea);
        std::cout<<"ea=exp(a):"<<std::endl;
        ch.print_matrix(ea,7);
        std::cout<<"niter: "<<niter<<std::endl;
        std::cout<<std::endl;

        //clean up memory used for cayley_hamilton tests:
        ch.delete_matrix(a);
        ch.delete_matrix(ea);
    }

    if(0) {
        std::cout<<"test #3:"<<std::endl;
        int n=2; // matrix size

        //initialize an instance of the sun_algebra class with matrix size n (i.e. Lie algebra su(n))
        chexp<ftype> ch(n);

        //allocate some matrices using utility functions from the cayley_hamilton sub-class of sun_algebra:
        ftype** a;
        ch.new_matrix(a,0);
        ftype** ea;
        ch.new_matrix(ea);

        a[0][0]=0.0;
        a[0][1]=1.0;
        a[1][0]=-39.0;
        a[1][1]=-40.0;

        //a[0][0]+=39.0;
        //a[1][1]+=39.0;


        // print out the matrix a[][]:
        std::cout<<"a:"<<std::endl;
        ch.print_matrix(a,7);

        // determine the matrix exponential ea[][] of a[][]:
        int niter=ch(a,ea);
        std::cout<<"ea=exp(a):"<<std::endl;
        //ch.matrix_mult_scalar(ea,std::exp(-39.0),a);
        //ch.print_matrix(a,7);
        ch.print_matrix(ea,7);
        std::cout<<"niter: "<<niter<<std::endl;
        std::cout<<std::endl;

        //clean up memory used for cayley_hamilton tests:
        ch.delete_matrix(a);
        ch.delete_matrix(ea);
    }

    if(0) {
        std::cout<<"test #4:"<<std::endl;
        int n=2; // matrix size

        //initialize an instance of the sun_algebra class with matrix size n (i.e. Lie algebra su(n))
        chexp<ftype> ch(n);

        //allocate some matrices using utility functions from the cayley_hamilton sub-class of sun_algebra:
        ftype** a;
        ch.new_matrix(a,0);
        ftype** ea;
        ch.new_matrix(ea);

        a[0][0]=-49.0;
        a[0][1]=24.0;
        a[1][0]=-64.0;
        a[1][1]=31.0;


        // print out the matrix a[][]:
        std::cout<<"a:"<<std::endl;
        ch.print_matrix(a,7);

        // determine the matrix exponential ea[][] of a[][]:
        int niter=ch(a,ea);
        std::cout<<"ea=exp(a):"<<std::endl;
        ch.print_matrix(ea,7);
        std::cout<<"niter: "<<niter<<std::endl;
        std::cout<<std::endl;

        //clean up memory used for cayley_hamilton tests:
        ch.delete_matrix(a);
        ch.delete_matrix(ea);
    }

    if(0) {
        std::cout<<"test #5:"<<std::endl;
        int n=4; // matrix size

        //initialize an instance of the sun_algebra class with matrix size n (i.e. Lie algebra su(n))
        chexp<ftype> ch(n);

        //allocate some matrices using utility functions from the cayley_hamilton sub-class of sun_algebra:
        ftype** a;
        ch.new_matrix(a,0);
        ftype** ea;
        ch.new_matrix(ea);

        a[0][1]=6.0;
        a[1][2]=6.0;
        a[2][3]=6.0;


        // print out the matrix a[][]:
        std::cout<<"a:"<<std::endl;
        ch.print_matrix(a,7);

        // determine the matrix exponential ea[][] of a[][]:
        int niter=ch(a,ea);
        std::cout<<"ea=exp(a):"<<std::endl;
        ch.print_matrix(ea,7);
        std::cout<<"niter: "<<niter<<std::endl;
        std::cout<<std::endl;

        //clean up memory used for cayley_hamilton tests:
        ch.delete_matrix(a);
        ch.delete_matrix(ea);
    }

    if(0) {
        std::cout<<"test #6:"<<std::endl;
        int n=2; // matrix size

        //initialize an instance of the sun_algebra class with matrix size n (i.e. Lie algebra su(n))
        chexp<ftype> ch(n);

        //allocate some matrices using utility functions from the cayley_hamilton sub-class of sun_algebra:
        ftype** a;
        ch.new_matrix(a,0);
        ftype** ea;
        ch.new_matrix(ea);

        a[0][0]=1.0;
        a[0][1]=1.0;
        a[1][1]=1.0;


        // print out the matrix a[][]:
        std::cout<<"a:"<<std::endl;
        ch.print_matrix(a,7);

        // determine the matrix exponential ea[][] of a[][]:
        int niter=ch(a,ea);
        std::cout<<"ea=exp(a):"<<std::endl;
        ch.print_matrix(ea,7);
        std::cout<<"niter: "<<niter<<std::endl;
        std::cout<<std::endl;

        //clean up memory used for cayley_hamilton tests:
        ch.delete_matrix(a);
        ch.delete_matrix(ea);
    }

    if(0) {
        std::cout<<"test #7:"<<std::endl;
        int n=2; // matrix size

        //initialize an instance of the sun_algebra class with matrix size n (i.e. Lie algebra su(n))
        chexp<ftype> ch(n);

        //allocate some matrices using utility functions from the cayley_hamilton sub-class of sun_algebra:
        ftype** a;
        ch.new_matrix(a,0);
        ftype** ea;
        ch.new_matrix(ea);

        a[0][0]=1.0+1.0e-5;
        a[0][1]=1.0;
        a[1][1]=1.0-1.0e-5;


        // print out the matrix a[][]:
        std::cout<<"a:"<<std::endl;
        ch.print_matrix(a,7);

        // determine the matrix exponential ea[][] of a[][]:
        int niter=ch(a,ea);
        std::cout<<"ea=exp(a):"<<std::endl;
        ch.print_matrix(ea,7);
        std::cout<<"niter: "<<niter<<std::endl;
        std::cout<<std::endl;

        //clean up memory used for cayley_hamilton tests:
        ch.delete_matrix(a);
        ch.delete_matrix(ea);
    }

    if(0) {
        std::cout<<"test #8:"<<std::endl;
        int n=3; // matrix size

        //initialize an instance of the sun_algebra class with matrix size n (i.e. Lie algebra su(n))
        chexp<ftype> ch(n);

        //allocate some matrices using utility functions from the cayley_hamilton sub-class of sun_algebra:
        ftype** a;
        ch.new_matrix(a,0);
        ftype** ea;
        ch.new_matrix(ea);

        a[0][0]=21.0;
        a[0][1]=17.0;
        a[0][2]=6.0;
        a[1][0]=-5.0;
        a[1][1]=-1.0;
        a[1][2]=-6.0;
        a[2][0]=4.0;
        a[2][1]=4.0;
        a[2][2]=16.0;

        // print out the matrix a[][]:
        std::cout<<"a:"<<std::endl;
        ch.print_matrix(a,7);

        // determine the matrix exponential ea[][] of a[][]:
        int niter=ch(a,ea);
        std::cout<<"ea=exp(a):"<<std::endl;
        ch.print_matrix_e(ea,7);
        std::cout<<"niter: "<<niter<<std::endl;
        std::cout<<std::endl;

        //clean up memory used for cayley_hamilton tests:
        ch.delete_matrix(a);
        ch.delete_matrix(ea);
    }

    if(0) {
        std::cout<<"test #9:"<<std::endl;
        int n=4; // matrix size

        //initialize an instance of the sun_algebra class with matrix size n (i.e. Lie algebra su(n))
        chexp<ftype> ch(n);

        //allocate some matrices using utility functions from the cayley_hamilton sub-class of sun_algebra:
        ftype** a;
        ch.new_matrix(a,0);
        ftype** ea;
        ch.new_matrix(ea);

        a[0][0]=1.0;
        a[0][1]=2.0;
        a[0][2]=2.0;
        a[0][3]=2.0;
        a[1][0]=3.0;
        a[1][1]=1.0;
        a[1][2]=1.0;
        a[1][3]=2.0;
        a[2][0]=3.0;
        a[2][1]=2.0;
        a[2][2]=1.0;
        a[2][3]=2.0;
        a[3][0]=3.0;
        a[3][1]=3.0;
        a[3][2]=3.0;
        a[3][3]=1.0;

        // print out the matrix a[][]:
        std::cout<<"a:"<<std::endl;
        ch.print_matrix(a,7);

        // determine the matrix exponential ea[][] of a[][]:
        int niter=ch(a,ea);
        std::cout<<"ea=exp(a):"<<std::endl;
        ch.print_matrix(ea,7);
        std::cout<<"niter: "<<niter<<std::endl;
        std::cout<<std::endl;

        //clean up memory used for cayley_hamilton tests:
        ch.delete_matrix(a);
        ch.delete_matrix(ea);
    }

    if(0) {
        std::cout<<"test #10:"<<std::endl;
        int n=3; // matrix size

        //initialize an instance of the sun_algebra class with matrix size n (i.e. Lie algebra su(n))
        chexp<ftype> ch(n);

        //allocate some matrices using utility functions from the cayley_hamilton sub-class of sun_algebra:
        ftype** a;
        ch.new_matrix(a,0);
        ftype** ea;
        ch.new_matrix(ea);

        a[0][0]=4.0;
        a[0][1]=2.0;
        a[0][2]=0.0;
        a[1][0]=1.0;
        a[1][1]=4.0;
        a[1][2]=1.0;
        a[2][0]=1.0;
        a[2][1]=1.0;
        a[2][2]=4.0;

        // print out the matrix a[][]:
        std::cout<<"a:"<<std::endl;
        ch.print_matrix(a,7);

        // determine the matrix exponential ea[][] of a[][]:
        int niter=ch(a,ea);
        std::cout<<"ea=exp(a):"<<std::endl;
        ch.print_matrix(ea,7);
        std::cout<<"niter: "<<niter<<std::endl;
        std::cout<<std::endl;

        //clean up memory used for cayley_hamilton tests:
        ch.delete_matrix(a);
        ch.delete_matrix(ea);
    }

    if(0) {
        std::cout<<"test #11:"<<std::endl;
        int n=3; // matrix size

        //initialize an instance of the sun_algebra class with matrix size n (i.e. Lie algebra su(n))
        chexp<ftype> ch(n);

        //allocate some matrices using utility functions from the cayley_hamilton sub-class of sun_algebra:
        ftype** a;
        ch.new_matrix(a,0);
        ftype** ea;
        ch.new_matrix(ea);

        a[0][0]=29.87942128909879;
        a[0][1]=0.7815750847907159;
        a[0][2]=-2.289519314033932;
        a[1][0]=0.7815750847907159;
        a[1][1]=25.72656945571064;
        a[1][2]=8.680737820540137;
        a[2][0]=-2.289519314033932;
        a[2][1]=8.680737820540137;
        a[2][2]=34.39400925519054;

        // print out the matrix a[][]:
        std::cout<<"a:"<<std::endl;
        ch.print_matrix(a,7);

        // determine the matrix exponential ea[][] of a[][]:
        int niter=ch(a,ea);
        std::cout<<"ea=exp(a):"<<std::endl;
        ch.print_matrix_e(ea,7);
        std::cout<<"niter: "<<niter<<std::endl;
        std::cout<<std::endl;

        //clean up memory used for cayley_hamilton tests:
        ch.delete_matrix(a);
        ch.delete_matrix(ea);
    }

    if(0) {
        std::cout<<"test #12:"<<std::endl;
        int n=3; // matrix size

        //initialize an instance of the sun_algebra class with matrix size n (i.e. Lie algebra su(n))
        chexp<ftype> ch(n);

        //allocate some matrices using utility functions from the cayley_hamilton sub-class of sun_algebra:
        ftype** a;
        ch.new_matrix(a,0);
        ftype** ea;
        ch.new_matrix(ea);

        a[0][0]=-131.0;
        a[0][1]=19.0;
        a[0][2]=18.0;
        a[1][0]=-390.0;
        a[1][1]=56.0;
        a[1][2]=54.0;
        a[2][0]=-387.0;
        a[2][1]=57.0;
        a[2][2]=52.0;


        // print out the matrix a[][]:
        std::cout<<"a:"<<std::endl;
        ch.print_matrix(a,7);

        // determine the matrix exponential ea[][] of a[][]:
        int niter=ch(a,ea);
        std::cout<<"ea=exp(a):"<<std::endl;
        ch.print_matrix(ea,7);
        std::cout<<"niter: "<<niter<<std::endl;
        std::cout<<std::endl;

        //clean up memory used for cayley_hamilton tests:
        ch.delete_matrix(a);
        ch.delete_matrix(ea);
    }

    if(0) {
        std::cout<<"test #13:"<<std::endl;
        int n=10; // matrix size

        //initialize an instance of the sun_algebra class with matrix size n (i.e. Lie algebra su(n))
        chexp<ftype> ch(n);

        //allocate some matrices using utility functions from the cayley_hamilton sub-class of sun_algebra:
        ftype** a;
        ch.new_matrix(a,0);
        ftype** ea;
        ch.new_matrix(ea);

        for(int i=0; i<n-1; ++i) {
            a[i][i+1]=1.0;
        }
        a[9][0]=1e-10;


        // print out the matrix a[][]:
        std::cout<<"a:"<<std::endl;
        ch.print_matrix_e(a,7);

        // determine the matrix exponential ea[][] of a[][]:
        int niter=ch(a,ea);
        std::cout<<"ea=exp(a):"<<std::endl;
        ch.print_matrix_e(ea,7);
        std::cout<<"niter: "<<niter<<std::endl;
        std::cout<<std::endl;

        //clean up memory used for cayley_hamilton tests:
        ch.delete_matrix(a);
        ch.delete_matrix(ea);
    }

    if(0) {
        std::cout<<"test #14:"<<std::endl;
        int n=3; // matrix size

        //initialize an instance of the sun_algebra class with matrix size n (i.e. Lie algebra su(n))
        chexp<ftype> ch(n);

        //allocate some matrices using utility functions from the cayley_hamilton sub-class of sun_algebra:
        ftype** a;
        ch.new_matrix(a,0);
        ftype** ea;
        ch.new_matrix(ea);

        a[0][0]=0.0;
        a[0][1]=1.0e-08;
        a[0][2]=0.0;
        a[1][0]=-2.0e+10-2.0e+08/3.0;
        a[1][1]=-3.0;
        a[1][2]=2.0e+10;
        a[2][0]=200.0e+00/3.0;
        a[2][1]=0.0;
        a[2][2]=-200.0e+00/3.0;


        // print out the matrix a[][]:
        std::cout<<"a:"<<std::endl;
        ch.print_matrix_e(a,7);

        // determine the matrix exponential ea[][] of a[][]:
        int niter=ch(a,ea);
        std::cout<<"ea=exp(a):"<<std::endl;
        ch.print_matrix_e(ea,7);
        std::cout<<"niter: "<<niter<<std::endl;
        std::cout<<std::endl;

        //clean up memory used for cayley_hamilton tests:
        ch.delete_matrix(a);
        ch.delete_matrix(ea);
    }

    if(1) {
        std::deque<std::string> dirnamel;
        dirnamel.push_back("test_matrix_files/sun_r0_1pi/");
        dirnamel.push_back("test_matrix_files/sun_r0_2pi/");
        dirnamel.push_back("test_matrix_files/sun_r0_5pi/");
        dirnamel.push_back("test_matrix_files/sun_r0_10pi/");
        dirnamel.push_back("test_matrix_files/sun_r0_npi/");
        //dirnamel.push_back("test_matrix_files/sun/");
        //dirnamel.push_back("test_matrix_files/a1.5/");
        //dirnamel.push_back("test_matrix_files/a5n/");
        ftype sclim=(ftype)1.0;
        std::string header="exp type     mat. size      av. rerr.       med. rerr.       max. rerr.    av. order    av. nb      time\n";
        if(1) {
            std::cout<<std::endl;
            std::cout<<"========= chexp() ========="<<std::endl;
            std::string methodename="chexp";
            std::cout<<"using rescaling target value sclim="<<sclim<<std::endl;
            for(auto dirnit=dirnamel.begin(); dirnit!=dirnamel.end(); ++dirnit) {
                std::string dirname=*dirnit;
                std::cout<<"processing file: "<<dirname<<std::endl;
                std::ostringstream oss;
                DIR* dir;
                struct dirent* ent;
                if((dir=opendir(dirname.c_str()))!=NULL) {
                    std::map<int,std::string> filelist;
                    while((ent=readdir(dir))!=NULL) {
                        if(ent->d_name[0]!='.') {
                            oss.str("");
                            oss<<dirname<<ent->d_name;
                            std::ifstream cfin(oss.str().c_str(),std::ifstream::in|std::ifstream::binary);
                            int n=0;
                            int nmat=0;
                            if(cfin.good()) {
                                cfin.read(reinterpret_cast<char*>(&(n)),sizeof(int));
                                cfin.read(reinterpret_cast<char*>(&(nmat)),sizeof(int));
                                if(n>0) {
                                    filelist[n]=oss.str();
                                }
                            }
                            cfin.close();
                        }
                    }
                    closedir(dir);
                    printf(header.c_str());
                    for(auto pair=filelist.begin(); pair!=filelist.end(); ++pair) {
                        std::ifstream cfin(pair->second.c_str(),std::ifstream::in|std::ifstream::binary);
                        int n=0;
                        int nmat=0;
                        if(cfin.good()) {
                            cfin.read(reinterpret_cast<char*>(&n),sizeof(int));
                            cfin.read(reinterpret_cast<char*>(&nmat),sizeof(int));
                        }

                        if(n<2) {
                            continue;
                        }

                        chexp<ctype> mexp(n,sclim);
                        ctype*** al;
                        mexp.new_matrix_array(al,nmat);
                        ctype*** earl;
                        mexp.new_matrix_array(earl,nmat);
                        ctype** ea;
                        mexp.new_matrix(ea);
                        ctype** dea;
                        mexp.new_matrix(dea);

                        ftype prec=_fprec*n;
                        ftype deanorm,eanorm;
                        ftype avrdeanorm=0.0,maxdeanorm=0.0,meddeanorm=0.0,avniter=0.0,avnb=0.0;
                        ftype* deanorml=new ftype[nmat];

                        int ic1,ic2,niter,nb;
                        std::cout<<FIXED_FLOAT(8);
                        std::complex<double> temp;
                        for(int imat=0; imat<nmat; ++imat) {
                            for(ic1=0; ic1<n; ++ic1) {
                                for(ic2=0; ic2<n; ++ic2) {
                                    cfin.read(reinterpret_cast<char*>(&(temp)),sizeof(std::complex<double>));
                                    al[imat][ic1][ic2]=(ctype)temp;
                                }
                            }
                            for(ic1=0; ic1<n; ++ic1) {
                                for(ic2=0; ic2<n; ++ic2) {
                                    cfin.read(reinterpret_cast<char*>(&(temp)),sizeof(std::complex<double>));
                                    earl[imat][ic1][ic2]=(ctype)temp;
                                }
                            }
                            niter=mexp(al[imat],ea,&nb);
                            mexp.matrix_sub(ea,earl[imat],dea);

                            mexp.get_frobenius_norm(dea,deanorm);
                            mexp.get_frobenius_norm(earl[imat],eanorm);
                            avrdeanorm+=deanorm/eanorm;
                            deanorml[imat]=deanorm/eanorm;
                            avniter+=(ftype)niter;
                            avnb+=(ftype)nb;
                        }
                        cfin.close();

                        avrdeanorm/=nmat;
                        avniter/=nmat;
                        avnb/=nmat;
                        std::stable_sort(deanorml,deanorml+nmat);
                        meddeanorm=deanorml[nmat/2];
                        maxdeanorm=deanorml[nmat-1];

                        int nts=getmcount();
                        for(int i=0; i<100; ++i) {
                            for(int imat=0; imat<nmat; ++imat) {
                                niter=mexp(al[imat],ea,&nb);
                            }
                        }
                        int nte=getmspan(nts);
                        printf("%-12s    % 3d       %0.5e      %0.5e      %0.5e      %6.3f     % 6.3f    %6.3fs\n",methodename.c_str(),n,avrdeanorm,meddeanorm,maxdeanorm,avniter,avnb,(ftype)nte/1000.0);

                        mexp.delete_matrix_array(al,nmat);
                        mexp.delete_matrix_array(earl,nmat);
                        mexp.delete_matrix(ea);
                        mexp.delete_matrix(dea);
                        delete[] deanorml;
                    }

                }
            }
        }

        if(1) {
            std::cout<<std::endl;
            std::cout<<"========= ochexp() ========="<<std::endl;
            std::string methodename="ochexp";
            std::cout<<"using rescaling target value sclim="<<sclim<<std::endl;
            for(auto dirnit=dirnamel.begin(); dirnit!=dirnamel.end(); ++dirnit) {
                std::string dirname=*dirnit;
                std::cout<<"processing file: "<<dirname<<std::endl;
                std::ostringstream oss;
                DIR* dir;
                struct dirent* ent;
                if((dir=opendir(dirname.c_str()))!=NULL) {
                    std::map<int,std::string> filelist;
                    while((ent=readdir(dir))!=NULL) {
                        if(ent->d_name[0]!='.') {
                            oss.str("");
                            oss<<dirname<<ent->d_name;
                            std::ifstream cfin(oss.str().c_str(),std::ifstream::in|std::ifstream::binary);
                            int n=0;
                            int nmat=0;
                            if(cfin.good()) {
                                cfin.read(reinterpret_cast<char*>(&(n)),sizeof(int));
                                cfin.read(reinterpret_cast<char*>(&(nmat)),sizeof(int));
                                if(n>0) {
                                    filelist[n]=oss.str();
                                }
                            }
                            cfin.close();
                        }
                    }
                    closedir(dir);
                    printf(header.c_str());
                    for(auto pair=filelist.begin(); pair!=filelist.end(); ++pair) {
                        std::ifstream cfin(pair->second.c_str(),std::ifstream::in|std::ifstream::binary);
                        int n=0;
                        int nmat=0;
                        if(cfin.good()) {
                            cfin.read(reinterpret_cast<char*>(&n),sizeof(int));
                            cfin.read(reinterpret_cast<char*>(&nmat),sizeof(int));
                        }

                        if(n<2) {
                            continue;
                        }

                        cayley_hamilton<ctype> mexp(n);
                        ctype*** al;
                        mexp.new_matrix_array(al,nmat);
                        ctype*** earl;
                        mexp.new_matrix_array(earl,nmat);
                        ctype** ea;
                        mexp.new_matrix(ea);
                        ctype** dea;
                        mexp.new_matrix(dea);

                        ftype prec=_fprec*n;
                        ftype deanorm,eanorm;
                        ftype avrdeanorm=0.0,maxdeanorm=0.0,meddeanorm=0.0,avniter=0.0,avnb=0.0;
                        ftype* deanorml=new ftype[nmat];

                        int ic1,ic2,niter,nb;
                        std::cout<<FIXED_FLOAT(8);
                        std::complex<double> temp;
                        for(int imat=0; imat<nmat; ++imat) {
                            for(ic1=0; ic1<n; ++ic1) {
                                for(ic2=0; ic2<n; ++ic2) {
                                    cfin.read(reinterpret_cast<char*>(&(temp)),sizeof(std::complex<double>));
                                    al[imat][ic1][ic2]=(ctype)temp;
                                }
                            }
                            for(ic1=0; ic1<n; ++ic1) {
                                for(ic2=0; ic2<n; ++ic2) {
                                    cfin.read(reinterpret_cast<char*>(&(temp)),sizeof(std::complex<double>));
                                    earl[imat][ic1][ic2]=(ctype)temp;
                                }
                            }
                            ftype anorm,sfac=1.0;
                            mexp.get_frobenius_norm(al[imat],anorm);
                            nb=0;
                            while(anorm*sfac>=sclim) {
                                sfac*=0.5;
                                ++nb;
                            }
                            niter=mexp(al[imat],ea,sfac);
                            mexp.matrix_sub(ea,earl[imat],dea);

                            mexp.get_frobenius_norm(dea,deanorm);
                            mexp.get_frobenius_norm(earl[imat],eanorm);
                            avrdeanorm+=deanorm/eanorm;
                            deanorml[imat]=deanorm/eanorm;
                            avniter+=(ftype)niter;
                            avnb+=(ftype)nb;
                        }
                        cfin.close();

                        avrdeanorm/=nmat;
                        avniter/=nmat;
                        avnb/=nmat;
                        std::stable_sort(deanorml,deanorml+nmat);
                        meddeanorm=deanorml[nmat/2];
                        maxdeanorm=deanorml[nmat-1];


                        int nts=getmcount();
                        for(int i=0; i<100; ++i) {
                            for(int imat=0; imat<nmat; ++imat) {
                                ftype anorm,sfac=1.0;
                                nb=0;
                                mexp.get_frobenius_norm(al[imat],anorm);
                                while(anorm*sfac>=sclim) {
                                    sfac*=0.5;
                                    ++nb;
                                }
                                niter=mexp(al[imat],ea,sfac);
                            }
                        }
                        int nte=getmspan(nts);
                        printf("%-12s    % 3d       %0.5e      %0.5e      %0.5e      %6.3f     % 6.3f    %6.3fs\n",methodename.c_str(),n,avrdeanorm,meddeanorm,maxdeanorm,avniter,avnb,(ftype)nte/1000.0);

                        mexp.delete_matrix_array(al,nmat);
                        mexp.delete_matrix_array(earl,nmat);
                        mexp.delete_matrix(ea);
                        mexp.delete_matrix(dea);
                        delete[] deanorml;
                    }

                }
            }
        }

        if(1) {
            std::cout<<std::endl;
            std::cout<<"========= nvexp() ========="<<std::endl;
            std::string methodename="nvexp";
            std::cout<<"using rescaling target value sclim="<<sclim<<std::endl;
            for(auto dirnit=dirnamel.begin(); dirnit!=dirnamel.end(); ++dirnit) {
                std::string dirname=*dirnit;
                std::cout<<"processing file: "<<dirname<<std::endl;
                std::ostringstream oss;
                DIR* dir;
                struct dirent* ent;
                if((dir=opendir(dirname.c_str()))!=NULL) {
                    std::map<int,std::string> filelist;
                    while((ent=readdir(dir))!=NULL) {
                        if(ent->d_name[0]!='.') {
                            oss.str("");
                            oss<<dirname<<ent->d_name;
                            std::ifstream cfin(oss.str().c_str(),std::ifstream::in|std::ifstream::binary);
                            int n=0;
                            int nmat=0;
                            if(cfin.good()) {
                                cfin.read(reinterpret_cast<char*>(&(n)),sizeof(int));
                                cfin.read(reinterpret_cast<char*>(&(nmat)),sizeof(int));
                                if(n>0) {
                                    filelist[n]=oss.str();
                                }
                            }
                            cfin.close();
                        }
                    }
                    closedir(dir);
                    printf(header.c_str());
                    for(auto pair=filelist.begin(); pair!=filelist.end(); ++pair) {
                        std::ifstream cfin(pair->second.c_str(),std::ifstream::in|std::ifstream::binary);
                        int n=0;
                        int nmat=0;
                        if(cfin.good()) {
                            cfin.read(reinterpret_cast<char*>(&n),sizeof(int));
                            cfin.read(reinterpret_cast<char*>(&nmat),sizeof(int));
                        }

                        if(n<2) {
                            continue;
                        }

                        nvexp<ctype> mexp(n,sclim);
                        ctype*** al;
                        mexp.new_matrix_array(al,nmat);
                        ctype*** earl;
                        mexp.new_matrix_array(earl,nmat);
                        ctype** ea;
                        mexp.new_matrix(ea);
                        ctype** dea;
                        mexp.new_matrix(dea);

                        ftype prec=_fprec*n;
                        ftype deanorm,eanorm;

                        ftype avrdeanorm=0.0,maxdeanorm=0.0,meddeanorm=0.0,avniter=0.0,avnb=0.0;
                        ftype* deanorml=new ftype[nmat];

                        int ic1,ic2,niter,nb;
                        std::cout<<FIXED_FLOAT(8);
                        std::complex<double> temp;
                        for(int imat=0; imat<nmat; ++imat) {
                            for(ic1=0; ic1<n; ++ic1) {
                                for(ic2=0; ic2<n; ++ic2) {
                                    cfin.read(reinterpret_cast<char*>(&(temp)),sizeof(std::complex<double>));
                                    al[imat][ic1][ic2]=(ctype)temp;
                                }
                            }
                            for(ic1=0; ic1<n; ++ic1) {
                                for(ic2=0; ic2<n; ++ic2) {
                                    cfin.read(reinterpret_cast<char*>(&(temp)),sizeof(std::complex<double>));
                                    earl[imat][ic1][ic2]=(ctype)temp;
                                }
                            }
                            niter=mexp(al[imat],ea,&nb);
                            mexp.matrix_sub(ea,earl[imat],dea);

                            mexp.get_frobenius_norm(dea,deanorm);
                            mexp.get_frobenius_norm(ea,eanorm);
                            avrdeanorm+=deanorm/eanorm;
                            deanorml[imat]=deanorm/eanorm;
                            avniter+=(ftype)niter;
                            avnb+=(ftype)nb;
                        }
                        cfin.close();

                        avrdeanorm/=nmat;
                        avniter/=nmat;
                        avnb/=nmat;
                        std::stable_sort(deanorml,deanorml+nmat);
                        meddeanorm=deanorml[nmat/2];
                        maxdeanorm=deanorml[nmat-1];

                        int nts=getmcount();
                        for(int i=0; i<100; ++i) {
                            for(int imat=0; imat<nmat; ++imat) {
                                niter=mexp(al[imat],ea,&nb);
                            }
                        }
                        int nte=getmspan(nts);
                        printf("%-12s    % 3d       %0.5e      %0.5e      %0.5e      %6.3f     % 6.3f    %6.3fs\n",methodename.c_str(),n,avrdeanorm,meddeanorm,maxdeanorm,avniter,avnb,(ftype)nte/1000.0);

                        mexp.delete_matrix_array(al,nmat);
                        mexp.delete_matrix_array(earl,nmat);
                        mexp.delete_matrix(ea);
                        mexp.delete_matrix(dea);
                        delete[] deanorml;
                    }

                }
            }
        }
        std::cout<<std::endl<<std::endl;
    }

}
