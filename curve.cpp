#include "curve.h"

using namespace std;

//oi kampules
template <class T>
curve<T>::curve() {}

template <class T>
curve<T>::curve(string inp)
{
    vector<string> tokens;
    stringstream check1(inp);

    string intermediate;
    while (getline(check1, intermediate, '\t'))
    {
        tokens.push_back(intermediate); //tokens has token[0] = id, token[1] = megethos, token[2] = all points (space separated now)
    }
    id = tokens[0]; //to id einai to prwto
    //std::cout << "Id is " << id << '\n';
    num_of_points = atoi(tokens[1].c_str());
    //std::cout << "Points num is " << num_of_points << '\n';

    vector<string> coords;
    stringstream check2(tokens[2]);
    string intermediate2;
    // Tokenizing w.r.t. space ' '
    while (getline(check2, intermediate2, ' '))
    {
        coords.push_back(intermediate2); //coords[i] = "(x0," (i artio) & coords[j] = "y0)" (j peritto)
    }

    for (unsigned int i = 0; i < coords.size(); i += 2)
    {
        stringstream tool_x(coords[i]);
        stringstream tool_y(coords[i + 1]);
        string tmp;
        T coordx;
        T coordy;
        std::string::size_type sz; // alias of size_t
        //bool x_done = false;                                   //elegxei an exoume grapsei kai to x
        //bool y_done = false;                                   //elegxei an exoume grapsei kai to y
        curve_point<T> temp; //constructor xwris orismata
        //x
        tmp = tool_x.str().substr(1, tool_x.str().length() - 1); //x = -6.4227999999999996
        if (typeid(coordx) == typeid(double))
        {
            coordx = stod(tmp, &sz);
        }
        else if (typeid(coordx) == typeid(float))
        {
            coordx = stof(tmp, &sz);
        }
        else if (typeid(coordx) == typeid(int))
        {
            coordx = stoi(tmp);
        }
        //std::cout << coordx << '\n';
        temp.set_x(coordx);
        //std::cout << temp.get_x() << '\n';
        //y
        tmp = tool_y.str().substr(0, tool_y.str().length() - 1); //y = 53.288000000000004
        if (typeid(coordy) == typeid(double))
        {
            coordy = stod(tmp, &sz);
        }
        else if (typeid(coordy) == typeid(float))
        {
            coordy = stof(tmp, &sz);
        }
        else if (typeid(coordy) == typeid(int))
        {
            coordy = stoi(tmp);
        }
        ////std::cout << coordy << '\n';
        temp.set_y(coordy);
        //std::cout << temp.get_y() << '\n';
        //pair
        my_points.push_back(temp);
    }
    /*for(int i =0; i< my_points.size(); i++){
      std:: cout << "eimai h" << id << my_points[i].get_x() << "," << my_points[i].get_y() << "\n";

    }*/

}

template <class T>
curve<T>::curve(vector<curve_point<T>> cps, int pointsnum, string idd)
{
    my_points.clear();
    my_points.resize(pointsnum);

    num_of_points = pointsnum;
    id = idd;
    typename vector<curve_point<T>>::iterator it2 = cps.begin();
    for (typename vector<curve_point<T>>::iterator it = my_points.begin(); it != my_points.end(); ++it)
    {
        *it = *it2;
        ++it2;
    }
}

template <class T>
curve<T>::~curve<T>() {}

template <class T>
unsigned int curve<T>::get_size()
{
    //return num_of_points;
    return my_points.size();

}

template <class T>
unsigned int curve<T>::get_v_size()
{
    return my_points.size();
}

template <class T>
std::string curve<T>::get_id()
{
    return id;
}

template <class T>
unsigned int curve<T>::get_id_as_int()
{
    string tmp = id;
    //int my_id = 0;
    int my_id2 = atoi(tmp.c_str()); //oi kampyles exoun mono mia leksh(arithmo) ws id
    /*for (unsigned int i = 0; i < id.length(); i++)
    {
        if (isdigit(id[i]))
        {
            tmp = id.substr(i, id.length() - 1);
            my_id = atoi(tmp.c_str());
            break;
        }
    }*/
    return my_id2;
}

template <class T>
std::vector<curve_point<T>> curve<T>::get_points()
{
    return my_points;
}

template <class T>
void curve<T>::set_id(std::string s1) { id = s1; }

template <class T>
void curve<T>::set_points(std::vector<curve_point<T>> vc)
{
    my_points.clear();
    for (unsigned int i = 0; i < vc.size(); i++)
    {
        my_points.push_back(vc[i]);
    }

    /*typename vector<T>::iterator it2 = vc.begin();
    for (typename vector<T>::iterator it = my_points.begin(); it != my_points.end(); ++it)
    {
        *it = *it2;
        ++it2;
    }*/
    return;
}

template <class T>
void curve<T>::set_num_of_pnts(unsigned int nop)
{
    num_of_points = nop;
}

template class curve<float>;
template class curve<int>;
template class curve<double>;
