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
#include "cayley_hamilton.h"
#include "sun_algebra.h"

int main()
{
    // Run some tests on the sun_algebra and/or the cayley_hamilton class:
    std::random_device r;
    std::seed_seq seedsq{r(),r(),r(),r(),r(),r(),r(),r()};
    std::mt19937_64 rng(seedsq);
    rng.discard(10000-1);
    std::uniform_real_distribution<> dis(0,1.0);

    if(0) {
        //run tests on the sun_algebra class
        int n=15; // matrix size

        //initialize an instance of the sun_algebra class with matrix size n (i.e. Lie algebra su(n))
        sun_algebra sun(n);

        //allocate some matrices using utility functions from the cayley_hamilton sub-class of sun_algebra:
        ctype** a;
        sun.ch.new_matrix(a);
        ctype** ea;
        sun.ch.new_matrix(ea);

        // specify a random vector "v" of magnitude 2*alpha in space spanned by Lie algebra basis
        //ftype alpha=0.5; //value for which "v" points to point within injectivty domain of exponential map exp_{id} : su(n)=T_{id}SU(n) -> SU(n))
        ftype alpha=5.0; //value for which "v" points to point outside of injectivity domain of exponetial map exp_{id} : su(n)=T_{id}SU(n) -> SU(n))
        ftype* v=new ftype[sun.ngen];
        ftype vnorm=0;
        for(int i=0; i<sun.ngen; ++i) {
            v[i]=1.0-2.0*dis(rng);
            vnorm+=v[i]*v[i];
        }
        vnorm=std::sqrt(vnorm);
        for(int i=0; i<sun.ngen; ++i) {
            v[i]*=2.0*alpha/vnorm;
        }

        // get the anit-hermitian matrix representation of "v": a=i*v.T ,  where T={lambda_0 / 2, ..., lambda_ngen / 2}
        sun.get_alg_mat_ah(v,a);

        std::cout<<"a:"<<std::endl;
        for(int i=0; i<n; ++i) {
            for(int j=0; j<n; ++j) {
                std::cout<<a[i][j]<<" ";
            }
            std::cout<<std::endl;
        }
        std::cout<<std::endl;

        // get the group matrix ea=Exp(i*v.T) 
        sun.get_grp_mat(v,ea);

        std::cout<<"ea=exp(a):"<<std::endl;
        for(int i=0; i<n; ++i) {
            for(int j=0; j<n; ++j) {
                std::cout<<ea[i][j]<<" ";
            }
            std::cout<<std::endl;
        }
        std::cout<<std::endl;

        // take the log of ea by determining the vector v in Lie algebra space for which ea=Exp(i*v.T) 
        // (note that the result of this does not need to agree with the initial "v" if "alpha" has been chose sufficently large,
        // so that "v" points outside of the injectivty domain of the exponential map exp_{id} : su(n)=T_{id}SU(n) -> SU(n))
        sun.get_log(ea,v);

        // get the corresponding anti hermitian matrix a=i*v.T
        sun.get_alg_mat_ah(v,a);

        std::cout<<"lea=log(ea):"<<std::endl;
        for(int i=0; i<n; ++i) {
            for(int j=0; j<n; ++j) {
                std::cout<<a[i][j]<<" ";
            }
            std::cout<<std::endl;
        }
        std::cout<<std::endl;

        // get the group matrix ea=Exp(i*v.T)
        // (this should now agree again with the group matrix that was obtained with the original "v") 
        sun.get_grp_mat(v,ea);

        std::cout<<"ea=exp(lea):"<<std::endl;
        for(int i=0; i<n; ++i) {
            for(int j=0; j<n; ++j) {
                std::cout<<ea[i][j]<<" ";
            }
            std::cout<<std::endl;
        }
        std::cout<<std::endl;

        //clean up memory used for sun_algebra tests:
        sun.ch.delete_matrix(a);
        sun.ch.delete_matrix(ea);
        delete[] v;
    }

    if(1) {
        //run tests on the cayley_hamilton class
        int n=5; // matrix size

        //initialize an instance of the sun_algebra class with matrix size n (i.e. Lie algebra su(n))
        cayley_hamilton<ctype> ch(n);

        //allocate some matrices using utility functions from the cayley_hamilton sub-class of sun_algebra:
        ctype** a;
        ch.new_matrix(a);
        ctype** ea;
        ch.new_matrix(ea);

        //test cayley_hamilton with a specific matrix to compare with known results  
        a[0][0]=ctype(0,0.636104061890406358);
        a[0][1]=ctype(0.415109700587091623,0.362251550192636662);
        a[0][2]=ctype(0.262379162994905379,0.136152290260253661);
        a[0][3]=ctype(1.000000000000000000,0.302688009059414176);
        a[0][4]=ctype(0.220271715386512135,0.004099496443320860);
        a[1][0]=ctype(-0.415109700587091623,0.362251550192636662);
        a[1][1]=ctype(0,0.513824029910692648);
        a[1][2]=ctype(0.158541011002899968,0.212282712425434950);
        a[1][3]=ctype(0.207720347312592146,0.319752570609831199);
        a[1][4]=ctype(0.174920870086319735,0.373122397555067487);
        a[2][0]=ctype(-0.262379162994905379,0.136152290260253661);
        a[2][1]=ctype(-0.158541011002899968,0.212282712425434950);
        a[2][2]=ctype(0,-0.0573176593232183687);
        a[2][3]=ctype(0.455199379282777317,0.138027494037257226);
        a[2][4]=ctype(0.385696416370347461,0.182111326006295315);
        a[3][0]=ctype(-1.000000000000000000,0.302688009059414176);
        a[3][1]=ctype(-0.207720347312592146,0.319752570609831199);
        a[3][2]=ctype(-0.455199379282777317,0.138027494037257226);
        a[3][3]=ctype(0,-0.860541393898962115);
        a[3][4]=ctype(0.205832835511187251,0.238577182257813214);
        a[4][0]=ctype(-0.220271715386512135,0.004099496443320860);
        a[4][1]=ctype(-0.174920870086319735,0.373122397555067487);
        a[4][2]=ctype(-0.385696416370347461,0.182111326006295315);
        a[4][3]=ctype(-0.205832835511187251,0.238577182257813214);
        a[4][4]=ctype(0,-0.232069038578918335);

        std::cout<<"a:"<<std::endl;
        for(int i=0; i<n; ++i) {
            for(int j=0; j<n; ++j) {
                std::cout<<a[i][j]<<" ";
            }
            std::cout<<std::endl;
        }
        std::cout<<std::endl;

        // determine the matrix exponential ea[][] of a[][]:
        ch(a,ea);
        std::cout<<"ea=exp(a):"<<std::endl;
        for(int i=0; i<n; ++i) {
            for(int j=0; j<n; ++j) {
                std::cout<<ea[i][j]<<" ";
            }
            std::cout<<std::endl;
        }
        std::cout<<std::endl;

        // determine the derivative of the exponential ea[][] of a[][] in the direction of a[][] itself:
        ctype** tdea;
        ch.new_matrix(tdea);
        ch.matrix_copy(a,tdea);
        ch(a,ea,tdea);
        std::cout<<"tdea=dexp(a)(tdea):"<<std::endl;
        for(int i=0; i<n; ++i) {
            for(int j=0; j<n; ++j) {
                std::cout<<tdea[i][j]<<" ";
            }
            std::cout<<std::endl;
        }
        std::cout<<std::endl;
        ch.delete_matrix(tdea);



        // allocate memory for the differential matrices
        ctype**** dea=new ctype***[n];
        for(int ic1=0; ic1<n; ++ic1) {
            dea[ic1]=new ctype**[n];
            for(int ic2=0; ic2<n; ++ic2) {
                ch.new_matrix(dea[ic1][ic2]);
            }
        }

        // determine the matrix exponential ea[][] of a[][] and the corresponding set differentials dea[i][j][][]=dea[][]/da[i][j]:
        ch(a,ea,dea);

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
        ch(a,ea2,dea2,0.2);

        // print the results for the matrix exponentials:
        std::cout<<"ea=exp(a):"<<std::endl;
        for(int i=0; i<n; ++i) {
            for(int j=0; j<n; ++j) {
                std::cout<<ea[i][j]<<" ";
            }
            std::cout<<std::endl;
        }
        std::cout<<std::endl;

        std::cout<<"ea2=exp(a):"<<std::endl;
        for(int i=0; i<n; ++i) {
            for(int j=0; j<n; ++j) {
                std::cout<<ea2[i][j]<<" ";
            }
            std::cout<<std::endl;
        }
        std::cout<<std::endl;

        std::cout<<"ea-ea2:"<<std::endl;
        for(int i=0; i<n; ++i) {
            for(int j=0; j<n; ++j) {
                std::cout<<ea[i][j]-ea2[i][j]<<" ";
            }
            std::cout<<std::endl;
        }
        std::cout<<std::endl;

        // print the results for the derivative terms:
        for(int ic1=0; ic1<n; ++ic1) {
            for(int ic2=0; ic2<n; ++ic2) {
                std::cout<<"dea["<<ic1<<"]["<<ic2<<"]:"<<std::endl;
                for(int i=0; i<n; ++i) {
                    for(int j=0; j<n; ++j) {
                        std::cout<<"  "<<dea[ic1][ic2][i][j]<<" ";
                    }
                    std::cout<<std::endl;
                }
                std::cout<<std::endl;

                std::cout<<"dea2["<<ic1<<"]["<<ic2<<"]:"<<std::endl;
                for(int i=0; i<n; ++i) {
                    for(int j=0; j<n; ++j) {
                        std::cout<<"  "<<dea2[ic1][ic2][i][j]<<" ";
                    }
                    std::cout<<std::endl;
                }
                std::cout<<std::endl;

                std::cout<<"dea["<<ic1<<"]["<<ic2<<"]-dea2["<<ic1<<"]["<<ic2<<"]:"<<std::endl;
                for(int i=0; i<n; ++i) {
                    for(int j=0; j<n; ++j) {
                        std::cout<<"  "<<dea[ic1][ic2][i][j]-dea2[ic1][ic2][i][j]<<" ";
                    }
                    std::cout<<std::endl;
                }
                std::cout<<std::endl;
            }
        }

        //clean up memory used for cayley_hamilton tests:
        ch.delete_matrix(a);
        ch.delete_matrix(ea);

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


}
