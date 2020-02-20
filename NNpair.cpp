#include "NNpair.h"
#include <string>
#include <stdlib.h> //atoi
using namespace std;

NNpair::NNpair(std::string q, std::string p)
{
    setq_id(q);
    setp_id(p);
}

string NNpair::getq_id()
{
    return q_id;
}

int NNpair::getq_id_as_int()
{
    string tmp = q_id;
    int my_id = 0;
    for (unsigned int i = 0; i < q_id.length(); i++)
    {
        if (isdigit(q_id[i]))
        {
            tmp = q_id.substr(i, q_id.length() - 1);
            my_id = atoi(tmp.c_str());
        }
    }
    return my_id;
}

void NNpair::setq_id(std::string idd)
{
    q_id = idd;
}

string NNpair::getp_id()
{
    return p_id;
}

int NNpair::getp_id_as_int()
{
    string tmp = p_id;
    int my_id = 0;
    for (unsigned int i = 0; i < p_id.length(); i++)
    {
        if (isdigit(p_id[i]))
        {
            tmp = p_id.substr(i, p_id.length() - 1);
            my_id = atoi(tmp.c_str());
        }
    }
    return my_id;
}

void NNpair::setp_id(std::string idd)
{
    p_id = idd;
}

double NNpair::get_distance()
{
    return distance;
}

void NNpair::set_distance(double dis)
{
    distance = dis;
}
