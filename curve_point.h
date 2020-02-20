#ifndef _CURVE_POINT_H
#define _CURVE_POINT_H

#include <utility> //pair

//ta shmeia twn kampulwn
template <class T>
class curve_point
{
private:
    std::pair<T, T> point; //leitourgei mono gia 2D curves

public:
    curve_point<T>();
    curve_point<T>(T, T);
    ~curve_point<T>();
    T get_x();
    T get_y();
    void set_x(T);
    void set_y(T);
};

#endif