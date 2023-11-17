#include<iostream>
#include <random>
#include "cayley_hamilton.h"

int main()
{
    std::random_device r;
    std::seed_seq seedsq{r(),r(),r(),r(),r(),r(),r(),r()};
    std::mt19937_64 rng(seedsq);
    rng.discard(10000-1);
    std::uniform_real_distribution<> dis(0,1.0);

    int n=5;

    sun_algebra sun(n);

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

    // get the anit-hermitian matrix a=i*v.T ,  where T={lambda_0 / 2,...lambda_ngen / 2}
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

    // get the anti hermitian matrix a=i*v.T
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

    sun.ch.delete_matrix(a);
    sun.ch.delete_matrix(ea);
    delete[] v;
}
