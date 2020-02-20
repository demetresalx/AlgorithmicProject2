#ifndef _CURVE_HT_H
#define _CURVE_HT_H
#include "curve.h"
#include "g_funs.h"
#include "my_vector.h"
#include "utils.h"

template <class T>
class curve_ht
{
private:
  int size;
  std::vector<std::vector<std::pair<curve<T> *, long int>>> table; //pinakas (o hash table) apo pinakes (h lista ana bucket) apo zeugaria (zeugari dianusmatos kai timhs g)
  g_funs<T> my_g;

public:
  curve_ht<T>(int size_to_be, int k_to_be, int dimensions, double w_to_be);
  curve_ht<T>(){};
  ~curve_ht();
  void hash_vector(my_vector<T> *v, curve<T> *cu);
  std::vector<std::string> hash_query(my_vector<T> *q,  curve<T> *cu ,double radius, bool repetition);
  //int get_vector_bucket_number(my_vector v);
};

#endif
