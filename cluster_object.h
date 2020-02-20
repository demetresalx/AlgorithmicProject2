#ifndef _CLUSTER_OBJECT_
#define _CLUSTER_OBJECT_

#include "my_vector.h"
#include "curve.h"
#include <vector>
#include <unordered_map>

template <class T>
class cluster
{
private:
  my_vector<T> center; //tautizoume to cluster me to kentro tou, epilegoume olo to antikeimeno kai oxi dekth gt mporei na mhn yparxei sto dataset
  std::unordered_map<std::string, my_vector<T> * > set_of_points; //periexei deiktes pros ta shmeia. epilogh ths domhs map gia eukoloterh metatroph
public:
  cluster<T>(){}; //empty contructor, idanikos gia periptwseis poy den exoyme yparkto kentro
  ~cluster();
  cluster<T>(my_vector<T> * c); //vasei enos deikth se my_vector
  void set_center(my_vector<T> * c); //vasei enos deikth se my_vector
  bool center_is_real(); //self explanatory
  std::string get_center_id();
  std::vector<T> get_center_coords(); //epistrefei se morfh vector tis suntetagmenes kentrou
  my_vector<T> * get_center_ptr();
  void incorporate_point(my_vector<T> * p); //enswmatwnei ena shmeio sto cluster
  void discorporate_point(my_vector<T> * p); //apoxwrizetai ena shmeio
  std::unordered_map<std::string, my_vector<T> * > *get_set_of_points();
  void print_cluster();
};


//to idio alla gia kampyles
template <class T>
class curve_cluster
{
private:
  curve<T> center; //tautizoume to cluster me to kentro tou, epilegoume olo to antikeimeno kai oxi dekth gt mporei na mhn yparxei sto dataset
  std::unordered_map<std::string, curve<T> * > set_of_curves; //periexei deiktes pros ta shmeia. epilogh ths domhs map gia eukoloterh metatroph
public:
  curve_cluster<T>(){}; //empty contructor, idanikos gia periptwseis poy den exoyme yparkto kentro
  ~curve_cluster();
  curve_cluster<T>(curve<T> * c); //vasei enos deikth se kampylh
  void set_center(curve<T> * c); //vasei enos deikth se kampylh
  //bool center_is_real(); //self explanatory
  std::string get_center_id();
  std::vector<curve_point<T>> get_center_points(); //epistrefei ta shmeia poy apoteloyn thn kampylh-kentro
  curve<T> * get_center_ptr();
  curve<T> get_center();
  void incorporate_point(curve<T> * p); //enswmatwnei ena shmeio sto cluster
  void discorporate_point(curve<T> * p); //apoxwrizetai ena shmeio
  std::unordered_map<std::string, curve<T> * > *get_set_of_curves();
  void print_cluster();
};


#endif
