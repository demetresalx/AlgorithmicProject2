#ifndef _H_FUNS_H_
#define _H_FUNS_H_

#include <vector>
#include "my_vector.h"
#include <cmath> //floor
#include <iostream>
#include <random> // uniform_real_distribution

/*Here we have a class of h functions which, concat'ed produce a G function which is associated with a HT - and we make L of them!*/

template <class T>
class h_funs
{
private:
    long long int m = 4294967291; //2 ^ 32 - 5 de xwraei se mikrotero!
    double w; //paronomastis - cell size gia ta ai
    long long int M; //M = 2^(32/k) de xwraei se mikrotero!
    std::vector<double> sis; //ta si pou xaraktirizoun thn h
    int dimensions; //edw d = 128 but we never know
public:
    h_funs<T>(int k, int dimens, double w_to_be);
    ~h_funs();
    long int individual_comp(long int ai, int expon);//dinei ton kathe typo poy xreiazetai gia to teliko h, expon o ek8eths tou m
    long int actual_h_function(my_vector<T> x); //long int takes at least 32 bits = 4 bytes
};

long int our_mod(long int a, long long int b); //a mod b
long int mod_pow(long int b, int e, long long int m);//https://en.wikipedia.org/wiki/Modular_exponentiation#Memory-efficient_method



#endif
