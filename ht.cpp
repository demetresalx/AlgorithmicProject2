#include "ht.h"
#include "utils.h"

//HT_CELL
template <class T>
ht_cell<T>::ht_cell(my_vector<T> vector_in_cell, long int g_val_to_be)
{
    my_v.set_v(vector_in_cell.get_v());
    my_v.set_id(vector_in_cell.get_id());
    g_value = g_val_to_be;
}

template <class T>
my_vector<T> ht_cell<T>::get_vector()
{
    return my_v;
}

template <class T>
long int ht_cell<T>::get_g_value()
{
    return g_value;
}

//?????????????????????????????
template <class T>
void ht_cell<T>::set_vector(my_vector<T> my_v_to_be)
{
    my_v.set_v(my_v_to_be.get_v());
    my_v.set_id(my_v_to_be.get_id());
    return;
}

template <class T>
void ht_cell<T>::set_g_value(long int my_g_val_to_be)
{
    g_value = my_g_val_to_be;
    return;
}

template <class T>
ht_cell<T>::~ht_cell() {}

//HT

template <class T>
ht<T>::ht(int size_to_be, int k_to_be, int dimensions, double w_to_be)
{
    //g_funs<T> my_g_to_be(k_to_be, dimensions, w_to_be);
    my_g = g_funs<T>(k_to_be, dimensions, w_to_be);
    size = size_to_be;
    //std::unordered_map<int, ht_cell<T> > table_to_be(size);
    table.clear();
    table.resize(size_to_be);
    for (int i = 0; i < size_to_be; i++)
    {
        table[i].clear();
    }
}

template <class T>
ht<T>::~ht() {}

template <class T>
void ht<T>::hash_vector(my_vector<T> *v)
{
    /*std::cout << "etsi hto o pinax prin hasharw to " << v->get_id_as_int() <<":";
    for(int i =0; i< table.size(); i++)
      for(int j =0; j< table[i].size(); j++)
        std::cout << table[i][j].first->get_id();
    std::cout << "\n";*/

    long int keyv = my_g.actual_g_function(*v); //ypologise thn timh ths g gia to vector v
    long int modded_keyv = our_mod(keyv, size); //kane thn timh auti mod table size
    //ht_cell<T> tmp_ht_cell(v,keyv); //krata se ena ht_cell to vector & to actual g value
    std::pair<my_vector<T> *, long int> thepair;
    thepair.first = v;
    thepair.second = keyv;
    //table[modded_keyv] = tmp_ht_cell; //vale auto to ht_cell sto katallilo bucket tou table
    table[modded_keyv].push_back(thepair);
    //std::cout << "my id is " << thepair.first->get_id() << "\n";
}

template <class T>
std::vector<std::string> ht<T>::hash_query(my_vector<T> *q, double radius, bool repetition) //emploutismenh me elegxo radius
{
    long int keyv = my_g.actual_g_function(*q); //ypologise thn timh ths g gia to vector v
    long int modded_keyv = our_mod(keyv, size); //kane thn timh auti mod table size
    std::vector<std::string> this_HT_potential_neighbs;
    this_HT_potential_neighbs.clear();
    for (unsigned int i = 0; i < table[modded_keyv].size(); i++)
    {
        //std::cout << table[modded_keyv][i].second << "\n"; //DEN PAIRNEI TIMES ID KAI VEC SWSTA EDW
        if(repetition == false){ //prwth fora, tsekare 0 < d <= r
          if( (table[modded_keyv][i].second == keyv) && (manhattan_distance(table[modded_keyv][i].first->get_v(), q->get_v() ) <= radius)  )                                    //prepei na exoyn to idio g epishs
            {this_HT_potential_neighbs.push_back(table[modded_keyv][i].first->get_id());  /*std::cout << "mojanga";*/}//valto sth lista pithanwn geitonwn
        }
        else{ //kanoyme epanalhpsh, ara koita apo r_prev kai meta
          if( (table[modded_keyv][i].second == keyv) && (manhattan_distance(table[modded_keyv][i].first->get_v(), q->get_v() ) <= radius) && (manhattan_distance(table[modded_keyv][i].first->get_v(), q->get_v() ) > radius/2)  )                                    //prepei na exoyn to idio g epishs
            {this_HT_potential_neighbs.push_back(table[modded_keyv][i].first->get_id()); /*std::cout << "mojanga";*/} //valto sth lista pithanwn geitonwn
        }

        //std::cout << "my id is " << table[modded_keyv][i].first->get_id_as_int() << "\n";
    }
    /*std::cout << "brhka tous: ";
    std::cout << "oi ypopshfioi geitonoi moy einai : ";
    for(unsigned int yod = 0; yod < this_HT_potential_neighbs.size(); yod++)
      std::cout << this_HT_potential_neighbs[yod];
    std::cout << "\n";*/
    return this_HT_potential_neighbs;
}

template class ht<float>;
template class ht<int>;
template class ht<double>;
