#include "my_vector.h"

using namespace std;

template <class T>
my_vector<T>::my_vector(std::vector<T> v_to_be, std::string id_to_be)
{
    id = id_to_be;
    typename vector<T>::iterator it2 = v_to_be.begin();
    for (typename vector<T>::iterator it = vec.begin(); it != vec.end(); ++it)
    {
        *it = *it2;
        ++it2;
    }
}

//to pneuma tou constructor einai na pairnei to input string - seira sto arxeio eisodou kai na kanei swsta tis arxikopoihseis
template <class T>
my_vector<T>::my_vector(std::string inp)
{
    /*vector<string> tokens;
    stringstream check1(inp);

    string intermediate;
    // Tokenizing w.r.t. space ' '
    while (getline(check1, intermediate, ' '))
    {
        tokens.push_back(intermediate);
    }

    string the_id = tokens[0]; //to id einai to prwto
    my_vector<T>::set_id(the_id);

    for (unsigned int i = 1; i < tokens.size(); i++)
    {
        stringstream tool(tokens[i]);
        int dimens_i = 0;
        tool >> dimens_i;
        vec.push_back(dimens_i); //vazei sto vector kathe stoixeio-diastash tou dianusmatos apo to input string
    }*/
    char comtool[inp.length()];
    strcpy(comtool, inp.c_str());
    char delimiter[]=" \t\r\n\v\f\0";
    //cmnd = strtok(comtool, delimiter);
    char *token = std::strtok(comtool, delimiter);
    std::string the_id(token);
    //std::cout << the_id;
    my_vector<T>::set_id(the_id);
    token = std::strtok(NULL, delimiter);
    while (token != NULL) { //oi diastaseis tou dianusmatos
        //std::cout << token << '\n';
        string gi(token);
        stringstream tool(gi);
        T dimens_i = 0;
        tool >> dimens_i;
        vec.push_back(dimens_i); //vazei sto vector kathe stoixeio-diastash tou dianusmatos apo to input string
        token = std::strtok(NULL, delimiter);
    }
}

template <class T>
void my_vector<T>::set_id(string idd)
{
    id = idd;
    return;
}

template <class T>
void my_vector<T>::set_v(std::vector<T> vv)
{
    vec.clear();
    for (unsigned int i = 0; i < vv.size(); i++)
    {
        vec.push_back(vv[i]);
    }

    /*    typename vector<T>::iterator it2 = vv.begin();
    for (typename vector<T>::iterator it = vec.begin(); it != vec.end(); ++it)
    {
        *it = *it2;
        ++it2;
    }*/
    if (vv.size() != vec.size())
    {
        fprintf(stderr, "Error in set_v\n");
        exit(-1);
    }

    return;
}

template <class T>
string my_vector<T>::get_id()
{
    return id;
}

template <class T>
int my_vector<T>::get_id_as_int()
{
    string tmp = id;
    int my_id = 0;
    for (unsigned int i = 0; i < id.length(); i++)
    {
        if (isdigit(id[i]))
        {
            tmp = id.substr(i, id.length() - 1);
            my_id = atoi(tmp.c_str());
            break;
        }
    }
    return my_id;
    //gia to input ths ekfwnhshs, poy htan kai sthn 1h ergasia
    /*string tmp = id;
    int my_id2 = atoi(tmp.c_str()); //exoun mono mia leksh(arithmo) ws id
    return my_id2;*/
}

template <class T>
vector<T> my_vector<T>::get_v()
{
    return vec;
}

template <class T>
my_vector<T>::~my_vector() {}

template class my_vector<float>;
template class my_vector<int>;
template class my_vector<double>;
