#include "cluster_object.h"
using namespace std;

template <class T>
cluster<T>::cluster(my_vector<T> * c){
  set_of_points.clear();
  set_center(c);
}


template <class T>
bool cluster<T>::center_is_real()
{
    if(center.get_id().length() == 0) //den einai alh8ino dianusma tou dataset
      return false;
    else
      return true;
}


template <class T>
void cluster<T>::set_center(my_vector<T> * c) //vasei enos deikth se my_vector
{
    center.set_id(c->get_id()); //to kanei idio
    center.set_v(c->get_v()); //ana8etei ton vector
}


template <class T>
std::string cluster<T>::get_center_id()
{
    return center.get_id();
}


template <class T>
vector<T> cluster<T>::get_center_coords()
{
    return center.get_v();
}

template <class T>
void cluster<T>::incorporate_point(my_vector<T> * p) //vasei enos deikth se my_vector
{
    set_of_points[p->get_id()] = p; //ana8esh timhs tou deikth sth swsth 8esh
}


template <class T>
void cluster<T>::discorporate_point(my_vector<T> * p) //vasei enos deikth se my_vector
{
    set_of_points.erase(p->get_id());  //diagrafei entelws to sugkekrimeno entry apo to cluster
}

template <class T>
cluster<T>::~cluster() {}

template <class T>
std::unordered_map<std::string, my_vector<T> * > *cluster<T>::get_set_of_points()
{
    return &set_of_points;
}

template <class T>
my_vector<T> * cluster<T>::get_center_ptr()
{
    return &center;
}

template <class T>
void cluster<T>::print_cluster()
{
    //std::cout << x.first << x.second.get_id() << "\n";
    std::cout << "Eimai to cluster " << get_center_id() << "\n";
    for (auto x : set_of_points)
    {
      std::cout << x.second->get_id() << " ";
    }
    std::cout << std::endl;


}

//SUNARTHSEIS CURVE_CLUSTER

template <class T>
curve_cluster<T>::~curve_cluster() {}


template <class T>
curve_cluster<T>::curve_cluster(curve<T> * c){
  set_of_curves.clear();
  set_center(c);
}


template <class T>
void curve_cluster<T>::set_center(curve<T> * c) //vasei enos deikth se my_vector
{
    center.set_id(c->get_id()); //to kanei idio
    center.set_points(c->get_points()); //ana8etei ton vector apo ta shmeia ths kampylhs
}


template <class T>
std::string curve_cluster<T>::get_center_id()
{
    return center.get_id();
}


template <class T>
std::vector<curve_point<T>> curve_cluster<T>::get_center_points()
{
    return center.get_points();
}


template <class T>
void curve_cluster<T>::incorporate_point(curve<T> *  p) //vasei enos deikth se my_vector
{
    set_of_curves[p->get_id()] = p; //ana8esh timhs tou deikth sth swsth 8esh
}


template <class T>
void curve_cluster<T>::discorporate_point(curve<T> *  p) //vasei enos deikth se my_vector
{
    set_of_curves.erase(p->get_id());  //diagrafei entelws to sugkekrimeno entry apo to cluster
}


template <class T>
std::unordered_map<std::string, curve<T> * > * curve_cluster<T>::get_set_of_curves()
{
    return &set_of_curves;
}

template <class T>
curve<T> * curve_cluster<T>::get_center_ptr()
{
    return &center;
}

template <class T>
curve<T> curve_cluster<T>::get_center()
{
    return center;
}


template <class T>
void curve_cluster<T>::print_cluster()
{
    //std::cout << x.first << x.second.get_id() << "\n";
    std::cout << "Eimai to cluster " << get_center_id() << "\n";
    for (auto x : set_of_curves)
    {
      std::cout << x.second->get_id() << " ";
    }
    std::cout << std::endl;


}



template class cluster<float>;
template class cluster<int>;
template class cluster<double>;
template class curve_cluster<float>;
template class curve_cluster<int>;
template class curve_cluster<double>;
