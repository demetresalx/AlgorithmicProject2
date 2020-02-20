#ifndef _G_FUNS_H
#define _G_FUNS_H

#include "h_funs.h"
#include <bitset>
#include <string>

/*Our class about the G-functions (k*h_functions) which produce the key which will determine where in the hash table our vector will go*/

template <class T>
class g_funs
{
private:
    int k;
    std::vector<h_funs<T>> my_h_funs;
public:
    g_funs<T>(int k_to_be, int dimensions, double w_to_be);
    g_funs<T>(){};
    ~g_funs();
    long int actual_g_function(my_vector<T> x);
};

#endif
