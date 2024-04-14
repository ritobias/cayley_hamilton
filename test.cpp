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

    if(1) {
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

    if(0) {
        //run tests on the sun_algebra class
        int n=7; // matrix size

        //initialize an instance of the sun_algebra class with matrix size n (i.e. Lie algebra su(n))
        cayley_hamilton<ctype> ch(n);

        //allocate some matrices using utility functions from the cayley_hamilton sub-class of sun_algebra:
        ctype** a;
        ch.new_matrix(a);
        ctype** ea;
        ch.new_matrix(ea);

        //test cayley_hamilton with a specific matrix to compare with known results  
        a[0][0]=ctype(1.74677858508888673,-0.23389258528841400);
        a[0][1]=ctype(-0.529248493063859533,0.400113921462462444);
        a[0][2]=ctype(0.681655621360357761,1.041694909948886898);
        a[0][3]=ctype(-0.463667528249608987,0.058761255602210655);
        a[0][4]=ctype(-0.442191549866464963,-0.502564502058318782);
        a[0][5]=ctype(-0.163999326830295742,-0.226221240781968779);
        a[0][6]=ctype(0.595390473720808518,0.670353101983073103);
        a[1][0]=ctype(-0.700671700449518125,0.498161436285685889);
        a[1][1]=ctype(0.67070417498161491,-1.41878434681873717);
        a[1][2]=ctype(-0.078867251861058321,-0.348136335940323430);
        a[1][3]=ctype(0.269910622598868321,-0.517449733024581309);
        a[1][4]=ctype(1.23028825345915115,0.70313582977058831);
        a[1][5]=ctype(-0.437917170337274741,-0.076853515874987783);
        a[1][6]=ctype(-0.681519222256972682,-0.676044345321808398);
        a[2][0]=ctype(-0.375766980511573515,0.363106141867388010);
        a[2][1]=ctype(1.36201884593811964,-0.97783866954593039);
        a[2][2]=ctype(1.278448503520886045,0.312041496682990476);
        a[2][3]=ctype(-0.322584264791542111,0.385664794789477367);
        a[2][4]=ctype(0.781643089830021669,-0.811435427297967561);
        a[2][5]=ctype(-0.139456179837426843,-0.748553183703038408);
        a[2][6]=ctype(0.1208450403428390280,0.0461175160484342597);
        a[3][0]=ctype(0.877895617789393950,0.069497701767636436);
        a[3][1]=ctype(0.250847078253968691,0.164740353806830889);
        a[3][2]=ctype(-0.514473115316049871,0.484978189521951340);
        a[3][3]=ctype(1.65102437762298697,0.30284317616689173);
        a[3][4]=ctype(-0.780072890260709757,-0.807846505334177529);
        a[3][5]=ctype(0.713645781797253660,-0.304466900261165410);
        a[3][6]=ctype(-0.0034425681260104585,0.0614707661443832574);
        a[4][0]=ctype(1.31211655303233413,-0.95532062834753811);
        a[4][1]=ctype(0.435602511392415308,0.109858143514514767);
        a[4][2]=ctype(-1.095365117645347418,-0.099790953144865254);
        a[4][3]=ctype(-0.875136429695418297,-0.499539836979216141);
        a[4][4]=ctype(1.58368198771569398,-0.33072994792154932);
        a[4][5]=ctype(-0.298858435310117864,0.088214679409882888);
        a[4][6]=ctype(0.731808573844583908,-0.127439297021755295);
        a[5][0]=ctype(0.077009296969940544,0.383459783252251116);
        a[5][1]=ctype(0.184589661006479904,0.462923111439274218);
        a[5][2]=ctype(-0.282949045600631518,0.614553703566781532);
        a[5][3]=ctype(0.15543621863398815,-1.56058139964717206);
        a[5][4]=ctype(-0.821631970961340708,-0.396713787673362353);
        a[5][5]=ctype(0.921492661221847190,-0.802899301238654477);
        a[5][6]=ctype(0.041882539794203765,-0.622478798915832024);
        a[6][0]=ctype(0.261871830034518760,-0.218851490311288807);
        a[6][1]=ctype(-0.718312747227648479,0.072033957715335456);
        a[6][2]=ctype(0.104297145171168619,0.461543103128907653);
        a[6][3]=ctype(0.750837156529914551,-0.338927830736385881);
        a[6][4]=ctype(-0.057904042562934785,0.244031236446473581);
        a[6][5]=ctype(-1.78867685506892153,-0.20898255117164019);
        a[6][6]=ctype(1.64888140330886664,-0.43549632210744375);



        std::cout<<"a:"<<std::endl;
        for(int i=0; i<n; ++i) {
            for(int j=0; j<n; ++j) {
                std::cout<<a[i][j]<<" ";
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
        std::cout<<"ea=exp(lea):"<<std::endl;
        for(int i=0; i<n; ++i) {
            for(int j=0; j<n; ++j) {
                std::cout<<ea[i][j]<<" ";
            }
            std::cout<<std::endl;
        }
        std::cout<<std::endl;

        std::cout<<"ea2=exp(lea):"<<std::endl;
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
