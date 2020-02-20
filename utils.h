#ifndef _UTILS_H_
#define _UTILS_H_

#include <stdio.h>
#include <iostream>
#include <vector>
#include <string>
#include <typeinfo>
#include <stdlib.h> //labs, llabs
#include <sstream>
#include <cmath> //abs with overload for float
#include <limits>
#include "cluster_object.h"
#include "my_vector.h"
#include "curve_point.h"
#include "curve.h"

/*Manhattan Distance*/
template <typename T>
double manhattan_distance(std::vector<T> v1, std::vector<T> v2);

template <typename T> /*template <class//typename T>*/
double manhattan_distance(std::vector<T> v1, std::vector<T> v2)
{ //sum twn (apoluti timi twn diaforwn twn  (suntetagmeni i tou v1, syntetagmeni i tou v2) )
	if (v1.size() != v2.size())
	{
		std::cerr << "Can't calculate distance, need same dimensions! Aborting...\n";
		exit(-1);
	}
	T temp;
	if ((typeid(temp) == typeid(float)) || (typeid(temp) == typeid(double)))
	{
		temp = 0.0;
	}
	else if (typeid(temp) == typeid(int))
	{
		temp = 0;
	}

	std::vector<T> diffs;
	for (unsigned int i = 0; i < v1.size(); ++i)
	{
		temp = abs(v1[i] - v2[i]);
		diffs.push_back(temp);
		//ws edw leitourgei ok
	}
	double result = 0.0;

	for (typename std::vector<T>::iterator it = diffs.begin(); it != diffs.end(); ++it)
	{
		//std::cerr << *it << "\n";
		result = result + (double)*it; //to athroisma twn apol. diaforwn
									   //std::cerr << result << "\n";
	}
	//std::cerr << "Ara telika exw " << result << "\n";
	return result;
}


/*Functions for DTW*/

template <typename T>
double true_euclidean(curve_point<T> cp1, curve_point<T> cp2)
{
  T x1 = cp1.get_x();
  T y1 = cp1.get_y();

  T x2 = cp2.get_x();
  T y2 = cp2.get_y();

  T x = x1 - x2; //calculating number to square in next step
  T y = y1 - y2;
  double dist;

  if (typeid(x) != typeid(y))
  {
    std::cerr << "Incompatible curve point types (What are you doing?!)\n";
    exit(-1);
  }
  if (typeid(x) == typeid(int))
  {
    int x_faux = x;
    int y_faux = y;
    dist = pow((double)x_faux, 2) + pow((double)y_faux, 2); //calculating Euclidean distance
    dist = sqrt(dist);
  }
  else if (typeid(x) == typeid(double))
  {
    double x_faux = x;
    double y_faux = y;
    dist = pow(x_faux, 2) + pow(y_faux, 2); //calculating Euclidean distance
    dist = sqrt(dist);
  }

  return dist;
}

//dtw me dynamic programming
template <typename T>
double dtw(curve<T> *c1, curve<T> *c2){

  double dist_mat[c1->get_size()+1][c2->get_size()+1];
  for(unsigned int i=1; i<= c1->get_size(); i++)
    dist_mat[i][0] = std::numeric_limits<double>::max();
  for(unsigned int i=1; i<= c2->get_size(); i++)
    dist_mat[0][i] = std::numeric_limits<double>::max();
  dist_mat[0][0] = 0.0;
  for(unsigned int i=1; i<= c1->get_size(); i++){
    for(unsigned int j=1; j<= c2->get_size(); j++){
      double mintrash = std::min(dist_mat[i-1][j-1], dist_mat[i][j-1]);
      mintrash = std::min(mintrash, dist_mat[i-1][j]);
      dist_mat[i][j] = true_euclidean(c1->get_points()[i-1], c2->get_points()[j-1] ) + mintrash;
    }
  }
  return dist_mat[c1->get_size()][c2->get_size()];

}


//format clusters based on kentra pou brethikan parapanw
template <typename T>
void format_clusters (std::vector<my_vector<T>> *cluster_centers, std::vector<cluster<T>> *clusters) {
  (*clusters).clear();
  for(unsigned int i=0; i<(*cluster_centers).size(); i++){
      cluster<T> onecatatime(&((*cluster_centers)[i]));
      (*clusters).push_back(onecatatime);
    }
}

template <typename T>
void format_curve_clusters (std::vector<curve<T>> *cluster_centers, std::vector<curve_cluster<T>> *clusters) {
  (*clusters).clear();
  for(unsigned int i=0; i<(*cluster_centers).size(); i++){
      curve_cluster<T> onecatatime(&((*cluster_centers)[i]));
      (*clusters).push_back(onecatatime);
    }
}


//shmantiko gia lsh kampylwn
template <typename T>
std::pair<double, T> calculate_delta(std::unordered_map<std::string, curve<T> > *kurv_array){
	std::vector<curve<T>> curves_array;
	for(auto x:*kurv_array){
		curves_array.push_back(x.second);
	}

	double delta = 0.0;                                        //mesi apostasi shmeiwn kampulws
    double inf = std::numeric_limits<double>::max();           //apeiro kai kala
    double max_coord = -1 * inf;                               //arxika -apeiro
    ;                                                          //h megisti metabliti pou uparxei sto dataset
    for (unsigned int i = 0; i < curves_array.size(); i++) //an theloume ola & query, sbhnoume to "-q" kai ola popa
    {
      double tmp = 0.0;
      double plithos_athroismatwn = 0.0;
      std::vector<curve_point<T>> shmeia = curves_array[i].get_points(); //pairnw to vector apo shmeia kathe kampulhs
      if (shmeia.size() < 2)
      {
        if (shmeia[i].get_x() > max_coord)
        {
          max_coord = shmeia[i].get_x();
        }
        if (shmeia[i].get_y() > max_coord)
        {
          max_coord = shmeia[i].get_y();
        }
        continue;
      }
      else
      {
        if (shmeia[i].get_x() > max_coord)
        {
          max_coord = shmeia[i].get_x();
        }
        if (shmeia[i].get_y() > max_coord)
        {
          max_coord = shmeia[i].get_y();
        }
        for (unsigned int j = 0; j < shmeia.size() - 1; j++) //pairnw kathe shmeio tou vector ^
        {
          tmp += true_euclidean<double>(shmeia[j], shmeia[j + 1]);
          plithos_athroismatwn += 1.0;
        }
        tmp = tmp / plithos_athroismatwn;
        delta += tmp;
        //std::cerr << "Tmp delta is " << delta << '\n';
      }
    }
    delta = delta / (double)(curves_array.size()); //ki edw sbhnw to "-q" an eimai stin periptosi pou thelw kai ta query
		std::pair<double, T> twoelems;
		twoelems.first = delta;
		twoelems.second = max_coord;
		return twoelems;
}


//Silhouette - My_Vector
template <typename T>
std::vector<double> Silhouette(std::vector<cluster<T>>* clusters)
{

  std::vector<double> epistrepsima; //double == to avg si tou cluster autou


  std::vector<my_vector<T>> ta_kentra; //edw krataw ta kentra twn cluster mono gia na ta exw pio eukola parakatw
  my_vector<T> tmp_kentro;
  ta_kentra.clear();
  for (unsigned int i = 0; i < clusters->size(); i++)
  {
    tmp_kentro.set_id((*clusters)[i].get_center_id());
    tmp_kentro.set_v((*clusters)[i].get_center_coords());
    ta_kentra.push_back(tmp_kentro);
  }
  if (ta_kentra.size() != clusters->size())
  {
    std::cerr << "Silhouette error (1)\n";
    exit(-1);
  }

  std::vector<double> avg_si; //gia upologismo tou average ana cluster
  double si = 0.0;

  for (unsigned int i = 0; i < clusters->size(); i++)//gia kathe cluster
  {
    double min1 = std::numeric_limits<double>::max(); //apeiro
    avg_si.clear();
    si = 0.0;
    std::unordered_map<std::string, my_vector<T> * > * clust_points = (*clusters)[i].get_set_of_points();
    for (auto x:*clust_points)//gia kathe shmeio tou cluster
    {
      std::vector<double> apostaseis_i;
      apostaseis_i.clear();
      for (auto y:*clust_points)//ta alla shmeia tou cluster autou
      {
        if (x.first == y.first)
          continue;
        
        apostaseis_i.push_back(manhattan_distance(x.second->get_v(), y.second->get_v()));
      }
      double mesi_apostasi_ai;
      //avg dist apo shmeia idiou cluster
      for (unsigned int j = 0; j < apostaseis_i.size(); j++)
      {
        mesi_apostasi_ai += apostaseis_i[j];
      }
      mesi_apostasi_ai /= apostaseis_i.size(); //a(i)
      apostaseis_i.clear();

      //vres next closest cluster
      for (unsigned int j = 0; j < ta_kentra.size(); j++)
      {
        if ((ta_kentra[j].get_id() == (*clusters)[i].get_center_id())&& (ta_kentra[j].get_v() == (*clusters)[i].get_center_coords())) //an einai to idi kentro
          apostaseis_i.push_back(-1.0);
        else
          apostaseis_i.push_back(manhattan_distance(x.second->get_v(), ta_kentra[j].get_v()));
      }
      min1 = std::numeric_limits<double>::max(); //apeiro
      int index_next_best_clust_cnt = -1;
      for (unsigned int j = 0; j < apostaseis_i.size(); j++)
      {
        if (apostaseis_i[j]<0)
          continue;
        if (min1 > apostaseis_i[j])
        {
          min1 = apostaseis_i[j];
          index_next_best_clust_cnt = j;
        }
      }
      if (index_next_best_clust_cnt < 0)
      {
        std::cerr << "Error in Silhouette (2)" << '\n';
        exit(-1);
      }
      apostaseis_i.clear();
      double mesi_apostasi_bi;
      //to index_next_best_clust_cnt krataei thn thesi tou next closest cluster opote:
      //ta_kentra[index_next_best_clust_cnt] ==  to kentro pou mas noiazei, tr thelw na vrw se poio cluster afora
      for (unsigned int j = 0; j < clusters->size(); j++)
      {
        //sto katallilo cluster aka to next best
        if ((*clusters)[j].get_center_id() == ta_kentra[index_next_best_clust_cnt].get_id())
        {
          std::unordered_map<std::string, my_vector<T> * > * clust_points2 = (*clusters)[j].get_set_of_points();
          for (auto y:*clust_points2)//ta alla shmeia tou cluster autou
          {
            apostaseis_i.push_back(manhattan_distance(x.second->get_v(), y.second->get_v()));
          }
          //avg dist apo shmeia tou next best cluster
          for (unsigned int k = 0; k < apostaseis_i.size(); k++)
          {
            mesi_apostasi_bi += apostaseis_i[k];
          }
          mesi_apostasi_bi /= apostaseis_i.size(); //b(i)
        }
        else
        {
          continue;
        }

      }

      //s(i) = [b(i) - a(i)] / max{a(i), b(i)}
      si = mesi_apostasi_bi - mesi_apostasi_ai;
      if (mesi_apostasi_bi > mesi_apostasi_ai)
        si /= mesi_apostasi_bi;
      else
        si /= mesi_apostasi_ai;
      avg_si.push_back(si);
    }
    //avg s(i) per cluster
    si = 0.0;
    for (unsigned int j = 0; j < avg_si.size(); j++)
    {
      si += avg_si[j];
    }
    si /= avg_si.size();
    //std::cerr << si << '\n';
    epistrepsima.push_back(si);
  }
  //std::cerr << '\n';
  return epistrepsima;
}

//Silhouette - Curves
template <typename T>
std::vector<double> Silhouette_curve(std::vector<curve_cluster<T>>* clusters)
{

  std::vector<double> epistrepsima; //double == to avg si tou cluster autou


  std::vector<curve<T>> ta_kentra; //edw krataw ta kentra twn cluster mono gia na ta exw pio eukola parakatw
  curve<T> tmp_kentro;
  ta_kentra.clear();
  for (unsigned int i = 0; i < clusters->size(); i++)
  {
    tmp_kentro.set_id((*clusters)[i].get_center_id());
    tmp_kentro.set_points((*clusters)[i].get_center_points());
    ta_kentra.push_back(tmp_kentro);
  }
  if (ta_kentra.size() != clusters->size())
  {
    std::cerr << "Silhouette error (1)\n";
    exit(-1);
  }

  std::vector<double> avg_si; //gia upologismo tou average ana cluster
  double si = 0.0;

  for (unsigned int i = 0; i < clusters->size(); i++)//gia kathe cluster
  {
    double min1 = std::numeric_limits<double>::max(); //apeiro
    avg_si.clear();
    si = 0.0;
    std::unordered_map<std::string, curve<T> * > * clust_points = (*clusters)[i].get_set_of_curves();
    for (auto x:*clust_points)//gia kathe shmeio tou cluster
    {
      std::vector<double> apostaseis_i;
      apostaseis_i.clear();
      for (auto y:*clust_points)//ta alla shmeia tou cluster autou
      {
        if (x.first == y.first)
          continue;
        
        apostaseis_i.push_back(dtw(x.second, y.second));
      }
      double mesi_apostasi_ai;
      //avg dist apo shmeia idiou cluster
      for (unsigned int j = 0; j < apostaseis_i.size(); j++)
      {
        mesi_apostasi_ai += apostaseis_i[j];
      }
      mesi_apostasi_ai /= apostaseis_i.size(); //a(i)
      apostaseis_i.clear();

      //vres next closest cluster
      for (unsigned int j = 0; j < ta_kentra.size(); j++)
      {
        if (ta_kentra[j].get_id() == (*clusters)[i].get_center_id()) //an einai to idi kentro
          apostaseis_i.push_back(-1.0);
        else
          apostaseis_i.push_back(dtw(x.second, &ta_kentra[j]));
      }
      min1 = std::numeric_limits<double>::max(); //apeiro
      int index_next_best_clust_cnt = -1;
      for (unsigned int j = 0; j < apostaseis_i.size(); j++)
      {
        if (apostaseis_i[j]<0)
          continue;
        if (min1 > apostaseis_i[j])
        {
          min1 = apostaseis_i[j];
          index_next_best_clust_cnt = j;
        }
      }
      if (index_next_best_clust_cnt < 0)
      {
        std::cerr << "Error in Silhouette (2)" << '\n';
        exit(-1);
      }
      apostaseis_i.clear();
      double mesi_apostasi_bi;
      //to index_next_best_clust_cnt krataei thn thesi tou next closest cluster opote:
      //ta_kentra[index_next_best_clust_cnt] ==  to kentro pou mas noiazei, tr thelw na vrw se poio cluster afora
      for (unsigned int j = 0; j < clusters->size(); j++)
      {
        //sto katallilo cluster aka to next best
        if ((*clusters)[j].get_center_id() == ta_kentra[index_next_best_clust_cnt].get_id())
        {
          std::unordered_map<std::string, curve<T> * > * clust_points2 = (*clusters)[j].get_set_of_curves();
          for (auto y:*clust_points2)//ta alla shmeia tou cluster autou
          {
            apostaseis_i.push_back(dtw(x.second, y.second));
          }
          //avg dist apo shmeia tou next best cluster
          for (unsigned int k = 0; k < apostaseis_i.size(); k++)
          {
            mesi_apostasi_bi += apostaseis_i[k];
          }
          mesi_apostasi_bi /= apostaseis_i.size(); //b(i)
        }
        else
        {
          continue;
        }

      }

      //s(i) = [b(i) - a(i)] / max{a(i), b(i)}
      si = mesi_apostasi_bi - mesi_apostasi_ai;
      if (mesi_apostasi_bi > mesi_apostasi_ai)
        si /= mesi_apostasi_bi;
      else
        si /= mesi_apostasi_ai;
      avg_si.push_back(si);
    } //telos ana shmeio
    //avg s(i) per cluster
    si = 0.0;
    for (unsigned int j = 0; j < avg_si.size(); j++)
    {
      si += avg_si[j];
    }
    si /= avg_si.size();
    //std::cerr << si << '\n';
    epistrepsima.push_back(si);
  } //telos ana cluster
  //std::cerr << '\n';
  return epistrepsima;
}


//allo function gia oliko meso si (idio gia curves & vectors)
template <typename T>
double Silhouette_oliko(std::vector<double> sis)
{
  double finalf = 0.0;
  for (unsigned int i = 0; i < sis.size(); i++)
  {
    finalf += sis[i];
  }
  finalf /= sis.size();
  return finalf;
}
#endif
