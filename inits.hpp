#include "my_vector.h"
#include "curve.h"
#include <set>
#include <unordered_map>
#include <algorithm>
#include "cluster_object.h"
#include <random>
#include <cmath>

//random init
template <typename T>
void initialise_centers(int clusters, std::unordered_map<std::string, my_vector<T>> *vectors_array, std::vector<cluster<T>> * clusts) {
  //srand(time(NULL));
  std::set<int> ids; //indexes ston pinaka cluster
  //logw anwmalias do8entwn input files, prepei na kanoyme workaround auto to problhma...
  std::vector<std::string> keys;
  for(auto kv : *vectors_array)
    keys.push_back(kv.first);


  for (int i = 0; i < clusters; i++)
  { //epilegw tuxaio int anamesa sto 0 kai to n
    std::random_device rand_dev;
    std::mt19937 generator(rand_dev());
    std::uniform_real_distribution<double>  distr(0, keys.size());
    int rand_id = floor(distr(generator));
    //int rand_id = rand()%keys.size();
    while (ids.find(rand_id) != ids.end()) //uparxei hdh auto pou brhke h rand
    {
      rand_id = floor(distr(generator));
      //rand_id = rand()%keys.size();
    }
    ids.insert(rand_id); //en telei bazw to unique value sto set m
  }
  std::vector<int> v_ids(ids.begin(), ids.end()); //epistrefw to set alla se vector apo indexes
  std::vector<my_vector<T>> ta_kentra;
  ta_kentra.clear();
  for (int i = 0; i < clusters; i++)
  {
    my_vector<T> one_v_atime;
    one_v_atime.set_id((*vectors_array)[keys[v_ids[i]]].get_id());
    one_v_atime.set_v((*vectors_array)[keys[v_ids[i]]].get_v());
    ta_kentra.push_back(one_v_atime);
  }
  if (ta_kentra.size() != (unsigned int)clusters)
  {
    std::cerr << "Error in initialise_centers.\n";
    exit(-2);
  }
  format_clusters(&ta_kentra, clusts);
  //return ta_kentra;
}

//BINARY SEARCH gia to telos ths kmeans++
std::string binarySearch(std::vector<std::pair<std::string, double>> *arr, int l, int r, double x)
{
        int mid = l + (r - l) / 2;

        if((r - l) == 1)//kratame ton aristero se periptwsh poy teleiwsei anamesa se dyo
          return (*arr)[l].first;

        if((r-mid == 1) && (mid-l==1)){ //ean ta l,mid, r einai diadoxika elegxw autes tis 2 periptwseis
          if(x > (*arr)[mid].second)
            return (*arr)[mid].first;
          else
            return (*arr)[l].first;
        }

        if(x > (*arr)[r].second) //oriakh timh ektos pinaka
          return (*arr)[r].first;

        // If the element belongs to the middle
        // itself
        if ((*arr)[mid].second == x)
            return (*arr)[mid].first;

        // If element is smaller than mid, then
        // it can only belong to left subarray
        if ((*arr)[mid].second > x)
            return binarySearch(arr, l, mid - 1, x);

        // Else the element belongs
        // in right subarray
        return binarySearch(arr, mid + 1, r, x);

}

//kmeans++
template <typename T>
void initialise_centers_plus(int clusters, std::unordered_map<std::string, my_vector<T>> *vectors_array, std::vector<cluster<T>> * clusts){
  //logw anwmalias do8entwn input files, prepei na kanoyme workaround auto to problhma...
  std::vector<std::string> keys; //ta ids
  for(auto kv : *vectors_array)
    keys.push_back(kv.first);

  //distance matrix - bazw edw distance otan auti ipologistei
  std::unordered_map <std::string ,std::unordered_map<std::string, double>> Distance_Map;
  //double Distance_Matrix[keys.size()][keys.size()]; //all points with all points distances
  for (unsigned int i = 0; i < keys.size(); i++)
  {
    for (unsigned int j = 0; j < keys.size(); j++)
    {
      Distance_Map[keys[i]][keys[j]] =-1.0; //arxika -1 gt den gnt negative distance - 8a upologizw & 8a vazw ligo ligo o,ti m xreiazetai => less upologismoi
    }
  }

  std::vector<std::string> ids_kentrwn; //ta Ids twn kentrwn
  ids_kentrwn.clear();
  int metrhths = 0;

  //to 1o kentro random
  std::random_device rand_dev;
  std::mt19937 generator(rand_dev());
  std::uniform_real_distribution<double>  distr(0, keys.size());
  int rand_id_1 = floor(distr(generator));

  std::vector<my_vector<T>> ta_kentra; //pou 8a epistrepsoume
  ta_kentra.clear();

  ids_kentrwn.push_back(keys[rand_id_1]); //balame to 1o kentro
  metrhths++;

  //array P
  std::vector<std::pair<std::string, double>> P_array; //auto pou mporei na ginei BST later on
  P_array.clear();

  std::vector<double> apostaseis; //edw 8a valw tis apostaseis gia na vrw thn max kai na diairesw me auti tis distances pou 8a balw sto P_array
  apostaseis.clear();

  //twra ta upoloipa kentra
  while(metrhths < clusters){

    for (unsigned int j = 0; j < keys.size(); j++) //kathe my_vector
    {
      //elegxw an einai kentro idi
      if (std::find(ids_kentrwn.begin(), ids_kentrwn.end() , (*vectors_array)[keys[j]].get_id()) != ids_kentrwn.end()) //uparxei hdh auto to my_v ws kentro
        continue;

      //an den einai kentro:
      double min1 = std::numeric_limits<double>::max(); //apeiro
      for (unsigned int i = 0; i < ids_kentrwn.size(); i++) //vres thn apostasi tou j apo kathe kentro i
      {
        //elegxw an exw upologisei to dist panw
        //if ((Distance_Map[ids_kentrwn[i]][(*vectors_array)[keys[j]].get_id()] > 0) && (Distance_Map[(*vectors_array)[keys[j]].get_id()][ids_kentrwn[i]] > 0)) //exei upologistei
        if (Distance_Map[ids_kentrwn[i]][(*vectors_array)[keys[j]].get_id()] < 0) //den tin exw upologisei
        {
          Distance_Map[ids_kentrwn[i]][(*vectors_array)[keys[j]].get_id()] = manhattan_distance((*vectors_array)[ids_kentrwn[i]].get_v(), (*vectors_array)[keys[j]].get_v()); //upologizw to distance auto
        }

        if (Distance_Map[ids_kentrwn[i]][(*vectors_array)[keys[j]].get_id()] < min1)
        {
          min1 = Distance_Map[ids_kentrwn[i]][(*vectors_array)[keys[j]].get_id()]; //krata tin elaxisti apostasi
        }
      } //telos for gia kathe kentro poy exw mexri stigmhs
      //else einai ypologismeno hdh, pame na kanoume push back:
      std::pair<std::string, double> temp_pair;
      temp_pair.first = (*vectors_array)[keys[j]].get_id(); //se poio shmeio anaferomai?
      temp_pair.second = min1*min1;
      apostaseis.push_back(min1); //apothikeuw
      P_array.push_back(temp_pair);
    } //telos for gia kathe shmeio
    //find max in apostaseis
    double max1 = 0.0;
    for (unsigned int i = 0; i < apostaseis.size(); i++)
    {
      if (apostaseis[i] > max1)
      {
        max1 = apostaseis[i];
      }
    }

    //divide all pair.second me autin tin timi
    for (unsigned int i = 0; i < P_array.size(); i++)
    {
      P_array[i].second = P_array[i].second/max1; //normalized
    }

    for (unsigned int i = 0; i < P_array.size(); i++)
    {
      double sum = 0.0;
      for (unsigned int j = 0; j < i; j++)
      {
        sum += P_array[j].second;
      }
      P_array[i].second += sum;
    }
    //bres random x (float) anamesa sto 0 kai to P_array[last].second
    double megisto = P_array[P_array.size()-1].second; //last value is tops
    std::uniform_real_distribution<double>  distr2(0, megisto);
    double x = distr2(generator);
    //auto to x einai anamesa se alles 2 times P_array[i].second kai P_array[i+1].second
    //binary search - P_array is already sorted! :D
    std::string to_become_center;
    to_become_center = binarySearch(&P_array, 0, P_array.size()-1, x);
    ids_kentrwn.push_back(to_become_center);
    metrhths++;
    //push back to x sta kentra kai ksanamana
  }; //telos megalhs while poy dhmiourgei 1 kentro th fora

  for(unsigned int i = 0; i < ids_kentrwn.size(); i++){
    my_vector<T> one_v_atime;
    one_v_atime.set_id((*vectors_array)[ids_kentrwn[i]].get_id());
    one_v_atime.set_v((*vectors_array)[ids_kentrwn[i]].get_v());
    ta_kentra.push_back(one_v_atime);
  }
  format_clusters(&ta_kentra, clusts);
  //return ta_kentra;

}//telos sunarthshs




////////////CURVE FUNCTIONS//////////
//random init
template <typename T>
void initialise_centers_curve(int clusters, std::unordered_map<std::string, curve<T>> *curves_array, std::vector<curve_cluster<T>> * clusts) {

  std::set<int> ids; //indexes ston pinaka cluster
  //logw anwmalias do8entwn input files, prepei na kanoyme workaround auto to problhma...
  std::vector<std::string> keys;
  for(auto kv : *curves_array)
    keys.push_back(kv.first);


  for (int i = 0; i < clusters; i++)
  { //epilegw tuxaio int anamesa sto 0 kai to n
    std::random_device rand_dev;
    std::mt19937 generator(rand_dev());
    std::uniform_real_distribution<double>  distr(0, keys.size());
    int rand_id = floor(distr(generator));
    //int rand_id = rand()%keys.size();
    while (ids.find(rand_id) != ids.end()) //uparxei hdh auto pou brhke h rand
    {
      rand_id = floor(distr(generator));
      //rand_id = rand()%keys.size();
    }
    ids.insert(rand_id); //en telei bazw to unique value sto set m
  }
  std::vector<int> v_ids(ids.begin(), ids.end()); //epistrefw to set alla se vector apo indexes
  std::vector<curve<T>> ta_kentra;
  ta_kentra.clear();
  for (int i = 0; i < clusters; i++)
  {
    curve<T> one_v_atime;
    one_v_atime.set_id((*curves_array)[keys[v_ids[i]]].get_id());
    one_v_atime.set_points((*curves_array)[keys[v_ids[i]]].get_points());
    ta_kentra.push_back(one_v_atime);
  }
  if (ta_kentra.size() != (unsigned int)clusters)
  {
    std::cerr << "Error in initialise_centers curves.\n";
    exit(-2);
  }
  format_curve_clusters(&ta_kentra, clusts);

}




//kmeans++
template <typename T>
void initialise_centers_plus_curve(int clusters, std::unordered_map<std::string, curve<T>> *curves_array, std::vector<curve_cluster<T>> * clusts){
  //logw anwmalias do8entwn input files, prepei na kanoyme workaround auto to problhma...
  std::vector<std::string> keys; //ta ids
  for(auto kv : *curves_array)
    keys.push_back(kv.first);

  //distance matrix - bazw edw distance otan auti ipologistei
  std::unordered_map <std::string ,std::unordered_map<std::string, double>> Distance_Map;
  //double Distance_Matrix[keys.size()][keys.size()]; //all points with all points distances
  for (unsigned int i = 0; i < keys.size(); i++)
  {
    for (unsigned int j = 0; j < keys.size(); j++)
    {
      Distance_Map[keys[i]][keys[j]] =-1.0; //arxika -1 gt den gnt negative distance - 8a upologizw & 8a vazw ligo ligo o,ti m xreiazetai => less upologismoi
    }
  }

  std::vector<std::string> ids_kentrwn; //ta Ids twn kentrwn
  ids_kentrwn.clear();
  int metrhths = 0;

  //to 1o kentro random
  std::random_device rand_dev;
  std::mt19937 generator(rand_dev());
  std::uniform_real_distribution<double>  distr(0, keys.size());
  int rand_id_1 = floor(distr(generator));

  std::vector<curve<T>> ta_kentra; //pou 8a epistrepsoume
  ta_kentra.clear();

  ids_kentrwn.push_back(keys[rand_id_1]); //balame to 1o kentro
  metrhths++;

  //array P
  std::vector<std::pair<std::string, double>> P_array; //auto pou mporei na ginei BST later on
  P_array.clear();

  std::vector<double> apostaseis; //edw 8a valw tis apostaseis gia na vrw thn max kai na diairesw me auti tis distances pou 8a balw sto P_array
  apostaseis.clear();

  //twra ta upoloipa kentra
  while(metrhths < clusters){

    for (unsigned int j = 0; j < keys.size(); j++) //kathe my_vector
    {
      //elegxw an einai kentro idi
      if (std::find(ids_kentrwn.begin(), ids_kentrwn.end() , (*curves_array)[keys[j]].get_id()) != ids_kentrwn.end()) //uparxei hdh auto to my_v ws kentro
        continue;

      //an den einai kentro:
      double min1 = std::numeric_limits<double>::max(); //apeiro
      for (unsigned int i = 0; i < ids_kentrwn.size(); i++) //vres thn apostasi tou j apo kathe kentro i
      {
        //elegxw an exw upologisei to dist panw
        //if ((Distance_Map[ids_kentrwn[i]][(*vectors_array)[keys[j]].get_id()] > 0) && (Distance_Map[(*vectors_array)[keys[j]].get_id()][ids_kentrwn[i]] > 0)) //exei upologistei
        if (Distance_Map[ids_kentrwn[i]][(*curves_array)[keys[j]].get_id()] < 0) //den tin exw upologisei
        {
          Distance_Map[ids_kentrwn[i]][(*curves_array)[keys[j]].get_id()] = dtw(&((*curves_array)[ids_kentrwn[i]]), &((*curves_array)[keys[j]])); //upologizw to distance auto
        }

        if (Distance_Map[ids_kentrwn[i]][(*curves_array)[keys[j]].get_id()] < min1)
        {
          min1 = Distance_Map[ids_kentrwn[i]][(*curves_array)[keys[j]].get_id()]; //krata tin elaxisti apostasi
        }
      } //telos for gia kathe kentro poy exw mexri stigmhs
      //else einai ypologismeno hdh, pame na kanoume push back:
      std::pair<std::string, double> temp_pair;
      temp_pair.first = (*curves_array)[keys[j]].get_id(); //se poio shmeio anaferomai?
      temp_pair.second = min1*min1;
      apostaseis.push_back(min1); //apothikeuw
      P_array.push_back(temp_pair);
    } //telos for gia kathe shmeio
    //find max in apostaseis
    double max1 = 0.0;
    for (unsigned int i = 0; i < apostaseis.size(); i++)
    {
      if (apostaseis[i] > max1)
      {
        max1 = apostaseis[i];
      }
    }

    //divide all pair.second me autin tin timi
    for (unsigned int i = 0; i < P_array.size(); i++)
    {
      P_array[i].second = P_array[i].second/max1; //normalized
    }

    for (unsigned int i = 0; i < P_array.size(); i++)
    {
      double sum = 0.0;
      for (unsigned int j = 0; j < i; j++)
      {
        sum += P_array[j].second;
      }
      P_array[i].second += sum;
    }
    //bres random x (float) anamesa sto 0 kai to P_array[last].second
    double megisto = P_array[P_array.size()-1].second; //last value is tops
    std::uniform_real_distribution<double>  distr2(0, megisto);
    double x = distr2(generator);
    //auto to x einai anamesa se alles 2 times P_array[i].second kai P_array[i+1].second
    //binary search - P_array is already sorted! :D
    std::string to_become_center;
    to_become_center = binarySearch(&P_array, 0, P_array.size()-1, x);
    ids_kentrwn.push_back(to_become_center);
    metrhths++;
    //push back to x sta kentra kai ksanamana
  }; //telos megalhs while poy dhmiourgei 1 kentro th fora

  for(unsigned int i = 0; i < ids_kentrwn.size(); i++){
    curve<T> one_v_atime;
    one_v_atime.set_id((*curves_array)[ids_kentrwn[i]].get_id());
    one_v_atime.set_points((*curves_array)[ids_kentrwn[i]].get_points());
    ta_kentra.push_back(one_v_atime);
  }
  format_curve_clusters(&ta_kentra, clusts);
  //return ta_kentra;

}//telos sunarthshs