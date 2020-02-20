#include "grid.h"
//#include "my_vector.cpp"
template <class T>
void grid<T>::define_hash_table(int size_to_be, int k_to_be, int dimensions, int w_to_be){
  curve_ht<T> myht(size_to_be, k_to_be, dimensions, w_to_be);
  hash_table = myht;

}

template <class T>
grid<T>::grid(double d_to_be, int dimensions)
{
    delta = d_to_be;
    //uniform distr gia t
    std::random_device rand_dev;
    std::mt19937 generator(rand_dev());
    std::uniform_real_distribution<double> distr(0, dimensions); //distr(0, delta)
    for (int i = 0; i < dimensions; ++i)
    {
        t.push_back(distr(generator));
    }
}

template <class T>
grid<T>::~grid() {}

template <class T>
curve<double> grid<T>::gridify(curve<T> *c)
{
    curve<double> gc;
    gc.set_id(c->get_id()); //diathrw to id giati eimai magkas
    //gia kathe shmeio ths c briskw to a = floor(c[i])/delta, kanw a * delta gia na brw to coordinate (se kathe aksona)
    //std::set<curve_point<double>> uniques;
    std::set<std::pair<double, double>> uniques;
    for (unsigned int i = 0; i < c->get_size(); i++)
    {
        T x1 = c->get_points()[i].get_x();
        double x2 = floor((double)x1 / delta);
        x2 *= delta; //to x tou grid curve gia to shmeio c[i]
        x2 += t[0]; //h metatopish gia to x
        T y1 = c->get_points()[i].get_y();
        double y2 = floor((double)y1 / delta);
        y2 *= delta; //to y tou grid curve gia to shmeio c[i]
        y2 += t[1]; //h metatopish gia to y
        //curve_point<double> cp(x2, y2);
        std::pair<double, double> cp_rep;
        cp_rep.first = x2;
        cp_rep.second = y2;
        uniques.insert(cp_rep); //ara krataei monadika
    }
    std::vector<curve_point<double>> tmp;
    std::vector<std::pair<double, double>> temp2(uniques.begin(), uniques.end());
    for (unsigned int i = 0; i < temp2.size(); i++)
    {
        curve_point<double> cp(temp2[i].first, temp2[i].second);
        tmp.push_back(cp);
    }
    gc.set_points(tmp);
    gc.set_num_of_pnts(tmp.size());
    return gc;
}



template class grid<float>;
template class grid<int>;
template class grid<double>;
