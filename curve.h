#ifndef _CURVE_H_
#define _CURVE_H_

#include <string>
#include <sstream>
#include <vector>        //vectors
#include "curve_point.h" //their points
#include <stdlib.h>      //atoi
#include <cctype>        //isdigit
#include <cstddef>       // std::size_t
#include <iostream>
#include <typeinfo>
#include <cstring>

//oi kampules
template <class T>
class curve
{
private:
    std::vector<curve_point<T>> my_points;
    unsigned int num_of_points; //plithos shmeiwn twn kampulwn
    std::string id;

public:
    //std::vector<curve_point<T>> my_points;
    curve<T>();                                                            //default constructor
    curve<T>(std::string inp);                                             //apo string
    curve<T>(std::vector<curve_point<T>>, int pointsnum, std::string idd); //parametropoihmenos
    ~curve<T>();                                                           //destructor
    unsigned int get_size();                                               //epistrefei to posa shmeia exei
    unsigned int get_v_size();
    std::string get_id();                                                  //epistrefei to ID ws string
    unsigned int get_id_as_int();                                          //epistrefei to ID ws int
    std::vector<curve_point<T>> get_points();                              //epistrefei ta shmeia ths kampulhs se ena vector apo pairs
    void set_id(std::string);
    void set_points(std::vector<curve_point<T>>);
    void set_num_of_pnts(unsigned int);
};

#endif
