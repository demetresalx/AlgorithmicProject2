#include "curve_ht.h"

//HT

template <class T>
curve_ht<T>::curve_ht(int size_to_be, int k_to_be, int dimensions, double w_to_be){
    //g_funs<T> my_g_to_be(k_to_be, dimensions, w_to_be);
    my_g = g_funs<T>(k_to_be, dimensions, w_to_be);
    size = size_to_be;
    //std::unordered_map<int, ht_cell<T> > table_to_be(size);
    table.clear();
    table.resize(size_to_be);
    for (int i = 0; i < size_to_be; i++) {
        table[i].clear();
    }

}

template <class T>
curve_ht<T>::~curve_ht(){}

template <class T>
void curve_ht<T>::hash_vector(my_vector<T> *v, curve<T> *cu) {
    /*std::cout << "etsi hto o pinax prin hasharw to " << v->get_id_as_int() <<":";
    for(int i =0; i< table.size(); i++)
      for(int j =0; j< table[i].size(); j++)
        std::cout << table[i][j].first->get_id_as_int();
    std::cout << "\n";*/

    long int keyv = my_g.actual_g_function(*v); //ypologise thn timh ths g gia to vector v
    long int modded_keyv = our_mod(keyv, size); //kane thn timh auti mod table size
    //ht_cell<T> tmp_ht_cell(v,keyv); //krata se ena ht_cell to vector & to actual g value
    std::pair<curve<T> * , long int> thepair;
    thepair.first = cu;
    thepair.second = keyv;
    //table[modded_keyv] = tmp_ht_cell; //vale auto to ht_cell sto katallilo bucket tou table
    table[modded_keyv].push_back(thepair);
    //std::cout << "my id is " << thepair.first->get_id_as_int() << "\n";
}

template <class T>
std::vector<std::string> curve_ht<T>::hash_query(my_vector<T> *q, curve<T> *cu ,double radius, bool repetition) {
    long int keyv = my_g.actual_g_function(*q); //ypologise thn timh ths g gia to vector v
    long int modded_keyv = our_mod(keyv, size); //kane thn timh auti mod table size
    std::vector<std::string> this_HT_potential_neighbs;
    this_HT_potential_neighbs.clear();
    for(unsigned int i =0; i<table[modded_keyv].size(); i++){
      //std::cout << table[modded_keyv][i].first->get_id_as_int() << "\n";
      if(repetition == false){ //prwth fora, tsekare 0 < d <= r
        if( (table[modded_keyv][i].second == keyv) && (dtw(table[modded_keyv][i].first, cu ) <= radius)  )                                    //prepei na exoyn to idio g epishs
          {this_HT_potential_neighbs.push_back(table[modded_keyv][i].first->get_id());  /*std::cout << "mojanga";*/}//valto sth lista pithanwn geitonwn
      }
      else{ //kanoyme epanalhpsh, ara koita apo r_prev kai meta
        if( (table[modded_keyv][i].second == keyv) && (dtw(table[modded_keyv][i].first, cu ) <= radius) && (dtw(table[modded_keyv][i].first, cu ) > radius/2)  )                                    //prepei na exoyn to idio g epishs
          {this_HT_potential_neighbs.push_back(table[modded_keyv][i].first->get_id()); /*std::cout << "mojanga";*/} //valto sth lista pithanwn geitonwn
      }
    }
    /*std::cout << "brhka tous: ";
    std::cout << "oi ypopshfioi geitonoi moy einai : ";
    for(unsigned int yod = 0; yod < this_HT_potential_neighbs.size(); yod++)
      std::cout << this_HT_potential_neighbs[yod];
    std::cout << "\n";*/
    return this_HT_potential_neighbs;
}

template class curve_ht<double>;
template class curve_ht<int>;
template class curve_ht<float>;
