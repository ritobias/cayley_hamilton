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
        int n=7; // matrix size

        //initialize an instance of the sun_algebra class with matrix size n (i.e. Lie algebra su(n))
        cayley_hamilton<ctype> ch(n);

        //allocate some matrices using utility functions from the cayley_hamilton sub-class of sun_algebra:
        ctype** a;
        ch.new_matrix(a);
        ctype** ea;
        ch.new_matrix(ea);

        //test cayley_hamilton with a specific matrix to compare with known results  
        a[0][0]=ctype(2.23107606336097330,0.33556431388627588);
        a[0][1]=ctype(-0.174074182131794689,0.965677647076969895);
        a[0][2]=ctype(-0.11614874957033151,1.86276813530164750);
        a[0][3]=ctype(0.329202696055506170,0.012547201750045694);
        a[0][4]=ctype(-0.613180576156485136,-0.365793161540751010);
        a[0][5]=ctype(-0.485276968650446976,0.139876450340040543);
        a[0][6]=ctype(1.76677939476217549,-0.27316175108439165);
        a[1][0]=ctype(-0.478003177976257247,-0.161082897120776887);
        a[1][1]=ctype(0.842232235753108943,-0.761304609762550560);
        a[1][2]=ctype(0.385613962793068531,0.157997073359470441);
        a[1][3]=ctype(-1.275682385139676099,0.227744278890307846);
        a[1][4]=ctype(1.39765898164490365,-0.26003369720956917);
        a[1][5]=ctype(-0.252064730158434410,1.056710979618895884);
        a[1][6]=ctype(0.599918598364021818,0.264862493520818892);
        a[2][0]=ctype(0.216764409788068816,0.380658349452795893);
        a[2][1]=ctype(0.685369361779604309,-0.237775262576115962);
        a[2][2]=ctype(1.032722447830450976,-0.626122953514803078);
        a[2][3]=ctype(0.193409452107556779,-0.075082148902527265);
        a[2][4]=ctype(-0.0605859783686185894,0.0422007020411549546);
        a[2][5]=ctype(0.64813979296481071,-1.80701548215797478);
        a[2][6]=ctype(0.367938580799691561,0.899215895472617827);
        a[3][0]=ctype(-0.472424869165117201,-0.118035500755481163);
        a[3][1]=ctype(-0.246260841660216811,1.011844043971473892);
        a[3][2]=ctype(0.236120335751716528,-0.445379772646140512);
        a[3][3]=ctype(1.36052830683374472,1.67421969140732711);
        a[3][4]=ctype(0.117035643901947917,0.610271016973902108);
        a[3][5]=ctype(-0.138390770804171885,0.540928334980740047);
        a[3][6]=ctype(-0.361100887627488759,-1.268562462770851354);
        a[4][0]=ctype(0.716693830939111877,0.202767315148578494);
        a[4][1]=ctype(-0.585668501546576718,-0.642563747578296220);
        a[4][2]=ctype(0.334739247720530351,-0.153235977179082787);
        a[4][3]=ctype(-0.326348598664571667,0.548872610367843187);
        a[4][4]=ctype(0.541447350051210242,-0.068987049119573390);
        a[4][5]=ctype(-0.81877223888807166,1.19643824885694323);
        a[4][6]=ctype(-0.52700080834327041,-1.35621716930745918);
        a[5][0]=ctype(-0.312138565670150618,-0.484721274529108603);
        a[5][1]=ctype(-0.234907845255175358,0.018860782001755700);
        a[5][2]=ctype(-0.718845987189393765,0.185271102398213612);
        a[5][3]=ctype(-0.107196171074459751,0.509018354969034849);
        a[5][4]=ctype(-0.687430044635832180,0.941286070507134108);
        a[5][5]=ctype(2.05766867504063615,1.03208124035395207);
        a[5][6]=ctype(0.1006747679859107433,0.0147583701119662839);
        a[6][0]=ctype(-1.275491491468996068,-0.580321308600470882);
        a[6][1]=ctype(-0.996819290198883243,0.320608676692837668);
        a[6][2]=ctype(0.128947147443858841,0.395565035562697534);
        a[6][3]=ctype(-0.426774799510843419,-0.405853382811553915);
        a[6][4]=ctype(0.470688962372544201,0.802857227769200283);
        a[6][5]=ctype(-0.068777099361351091,-0.124549920617004507);
        a[6][6]=ctype(1.56078915162169498,0.13973657154200797);

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

        // determine the matrix exponential ea2[][] of a[][] and the corresponding set differentials dea2[i][j][][]=dea2[][]/da[i][j]
        // using input matrix rescaling:
        ch(a,ea2,dea2,1);

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
