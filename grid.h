#ifndef _GRID_H_
#define _GRID_H_

#include <vector>
#include "curve.h"
#include <random> // uniform_real_distribution
#include "my_vector.h"
#include <set>
#include "curve_ht.h"
#include <algorithm>
//#include "my_vector.cpp"

template <class T>
class grid
{
private:
    std::vector<double> t;
    double delta;
    //curve_ht<T> hash_table;
    //pithanws deutero antikeimeno edw gia kybo

public:
    curve_ht<T> hash_table;
    grid<T>(double, int);
    ~grid();
    curve<double> gridify(curve<T> *);
    void define_hash_table(int size_to_be, int k_to_be, int dimensions, int w_to_be);
};

#endif
