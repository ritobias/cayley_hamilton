#include<iostream>
#include<random>
#include "sun_algebra.h"

int main()
{
    // Run some tests on the sun_algebra and cayley_hamilton classes:
    std::random_device r;
    std::seed_seq seedsq{r(),r(),r(),r(),r(),r(),r(),r()};
    std::mt19937_64 rng(seedsq);
    rng.discard(10000-1);
    std::uniform_real_distribution<> dis(0,1.0);

    int n=7; // matrix size

    //initialize an instance of the sun_algebra class with matrix size n (i.e. Lie algebra su(n))
    sun_algebra sun(n);

    //allocate some matrices using utility functions from the cayley_hamilton sub-class of sun_algebra:
    ctype** a;
    sun.ch.new_matrix(a);
    ctype** ea;
    sun.ch.new_matrix(ea);

    // specify a random vector "v" of length 2*alpha in Lie algebra space
    ftype alpha=5.0;
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

    // get the anit-hermitian matrix a=i*v.T ,  where T={lambda_0 / 2, ..., lambda_ngen / 2}
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
    sun.log_ah(ea,v);

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
    sun.get_grp_mat(v,ea);

    std::cout<<"ea=exp(lea):"<<std::endl;
    for(int i=0; i<n; ++i) {
        for(int j=0; j<n; ++j) {
            std::cout<<ea[i][j]<<" ";
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;

    if(0) {
        //test cayley_hamilton with a specific matrix to compare with known result  
        a[0][0]=ctype(0.939897479157750058,0.159364591527384149);
        a[0][1]=ctype(0.609542856276176626,-0.409090788157241629);
        a[0][2]=ctype(-0.039059897431674262,0.517193215519405116);
        a[0][3]=ctype(0.40725052326516726,-2.01846981108155516);
        a[0][4]=ctype(0.189293831063409556,-0.468510465124830266);
        a[0][5]=ctype(-0.364303188273199483,1.196075345647117086);
        a[0][6]=ctype(0.208483819847761022,-0.030197636391871884);
        a[1][0]=ctype(0.136085372465963523,0.204214878793515550);
        a[1][1]=ctype(1.82088253757827834,-0.69318655216260838);
        a[1][2]=ctype(0.399022830853585737,0.855576320616568304);
        a[1][3]=ctype(0.269381110148073252,0.335596971392543412);
        a[1][4]=ctype(1.324521170373122792,0.097719775321094394);
        a[1][5]=ctype(-0.245543994945725993,-0.992410947185864661);
        a[1][6]=ctype(-0.649313180460195042,-0.059181097125996425);
        a[2][0]=ctype(0.369790573575637641,-0.945517146780218789);
        a[2][1]=ctype(-0.811546011509055415,0.173352091039376230);
        a[2][2]=ctype(0.375919781491310153,-0.849214780092822886);
        a[2][3]=ctype(0.055651345004950875,0.739703997891538152);
        a[2][4]=ctype(1.161648742374914133,0.366298231248701811);
        a[2][5]=ctype(-0.550638339371502677,-0.456537891675555410);
        a[2][6]=ctype(0.407783835121606891,-0.539420960579541076);
        a[3][0]=ctype(0.856280517290871829,-0.523896017033833374);
        a[3][1]=ctype(0.249315672997534248,-0.015220738547211917);
        a[3][2]=ctype(0.100894015210061856,0.340356066446914124);
        a[3][3]=ctype(1.318137153107333681,-0.269881481911962792);
        a[3][4]=ctype(-0.276026436761462112,0.159212052668554386);
        a[3][5]=ctype(0.826383972857243170,0.244273280389912053);
        a[3][6]=ctype(0.329159165940204438,0.564315714237826739);
        a[4][0]=ctype(1.31279359176855429,0.79163098815521436);
        a[4][1]=ctype(-0.531740615339792813,0.164748070192097044);
        a[4][2]=ctype(0.184990984489392082,0.115482527117499458);
        a[4][3]=ctype(0.917862670069802480,0.111133145190458058);
        a[4][4]=ctype(1.54394803539122787,-0.96396619843393386);
        a[4][5]=ctype(0.246665289419386235,-0.696081518582821820);
        a[4][6]=ctype(-1.34207954183614791,0.96912030257780103);
        a[5][0]=ctype(0.480824399954413425,0.392561818066786137);
        a[5][1]=ctype(-0.57978937213263027,1.31145469358987082);
        a[5][2]=ctype(-0.226221396318719950,-0.061863930116264419);
        a[5][3]=ctype(-0.0348828144001546586,-0.0672002835720141840);
        a[5][4]=ctype(-1.266257924988785054,-0.120761198266586459);
        a[5][5]=ctype(2.68756401555185650,0.16744635193042443);
        a[5][6]=ctype(-0.195802332936151958,0.187006587994376163);
        a[6][0]=ctype(0.128063601261061158,0.074690275837103619);
        a[6][1]=ctype(-0.265471848659404633,-0.005268117652720439);
        a[6][2]=ctype(-0.317183260876006561,-1.232336323873950158);
        a[6][3]=ctype(-0.717352863593895747,-0.381845076386106669);
        a[6][4]=ctype(1.71968542947112091,-0.08166216734641273);
        a[6][5]=ctype(0.977767407219326747,-0.303859055216179626);
        a[6][6]=ctype(0.904955688097058369,0.230710556237842452);

        // allocate memory for the differential matrices
        ctype**** dea=new ctype***[n];
        for(int ic1=0; ic1<n; ++ic1) {
            dea[ic1]=new ctype**[n];
            for(int ic2=0; ic2<n; ++ic2) {
                sun.ch.new_matrix(dea[ic1][ic2]);
            }
        }

        // determine the matrix exponential ea[][] of a[][] and the corresponding set differentials dea[i][j][][]=dea[][]/da[i][j]:
        sun.ch(a,ea,dea);

        // allocate more memory to compare the results obtained with and wouthout scaling of the input matrix (should yield the same results):
        ctype** ea2;
        sun.ch.new_matrix(ea2);
        ctype**** dea2=new ctype***[n];
        for(int ic1=0; ic1<n; ++ic1) {
            dea2[ic1]=new ctype**[n];
            for(int ic2=0; ic2<n; ++ic2) {
                sun.ch.new_matrix(dea2[ic1][ic2]);
            }
        }

        // determine the matrix exponential ea2[][] of a[][] and the corresponding set differentials dea2[i][j][][]=dea2[][]/da[i][j]
        // using input matrix rescaling:
        sun.ch(a,ea2,dea2,1);

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

        //clean up extra memory used for cayley_hamilton tests:
        for(int ic1=0; ic1<n; ++ic1) {
            for(int ic2=0; ic2<n; ++ic2) {
                sun.ch.delete_matrix(dea[ic1][ic2]);
            }
            delete[] dea[ic1];
        }
        delete[] dea;

        sun.ch.delete_matrix(ea2);

        for(int ic1=0; ic1<n; ++ic1) {
            for(int ic2=0; ic2<n; ++ic2) {
                sun.ch.delete_matrix(dea2[ic1][ic2]);
            }
            delete[] dea2[ic1];
        }
        delete[] dea2;
    }

    //clean up memory used for sun_algebra tests:
    sun.ch.delete_matrix(a);
    sun.ch.delete_matrix(ea);
    delete[] v;
}
