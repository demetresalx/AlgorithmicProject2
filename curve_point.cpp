#include "curve_point.h"

using namespace std;

template <class T>
curve_point<T>::curve_point() {}

template <class T>
curve_point<T>::curve_point(T x, T y)
{
    point.first = x;
    point.second = y;
}

template <class T>
curve_point<T>::~curve_point<T>() {}

template <class T>
T curve_point<T>::get_x() { return point.first; }

template <class T>
T curve_point<T>::get_y() { return point.second; }

template <class T>
void curve_point<T>::set_x(T x_to_be) { point.first = x_to_be; }

template <class T>
void curve_point<T>::set_y(T y_to_be) { point.second = y_to_be; }

template class curve_point<double>;
template class curve_point<float>;
template class curve_point<int>;
