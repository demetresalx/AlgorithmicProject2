#include "g_funs.h"

template <class T>
g_funs<T>::g_funs(int k_to_be, int dimensions, double w_to_be) {
    k = k_to_be;
    for (int i = 0; i < k; i++) {
        h_funs<T> temp_h_fun(k_to_be, dimensions, w_to_be);
        my_h_funs.push_back(temp_h_fun);
    }
}

template <class T>
g_funs<T>::~g_funs(){}

template <class T>
long int g_funs<T>::actual_g_function(my_vector<T> x) {//epistrefei int >= 0


    //to apotelesma to kanw mod Table_Size (blepe lsh.cpp)
    //auto kathorize pou bazw to v sto Hash Table m
    //std::vector<T> the_v = x.get_v();

    unsigned long g_value1 = 0.0;
    long int temp;
    for(int i =0; i<k; i++){
        temp = my_h_funs[i].actual_h_function(x); //pairnw tis k h_funs gia ena vector v
        std::bitset<32> x(temp);
        x = x << i+1; //tis kanw concat (me left shift...
        g_value1 += x.to_ulong();
    }
    return g_value1;
}

template class g_funs<int>;
template class g_funs<float>;
template class g_funs<double>;
