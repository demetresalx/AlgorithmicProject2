#include "utils.h"
#include <unordered_map>
#include "cluster_object.h"
#include "my_vector.h"
#include <limits>

template <typename T>
void update_pam(std::vector<cluster<T>> * clusts){
  std::vector<my_vector<T>> nea_kentra;
  for(unsigned int i=0; i< (*clusts).size(); i++){ //gia ka8e cluster
    double min1 = std::numeric_limits<double>::max(); //min pairnei timh apeiro
    my_vector<T> min_dist_point;
    std::unordered_map<std::string, my_vector<T> * > * clust_points = (*clusts)[i].get_set_of_points(); //ta stoixeia tou cluster
    for(auto x:*clust_points) { //gia kathe stoixeio tou cluster
      double mesi_apostasi = 0.0;
      //sugkrine gia kathe shmeio != kentro
      if (x.first == (*clusts)[i].get_center_id())
        continue; //den psaxnw to idi kentro
      //else
      for(auto y:*clust_points) { //gia kathe stoixeio tou cluster
        //skip myself
        if (y.first == x.first)
          continue;
        //bres mean tou dist apo ta alla shmeia
        mesi_apostasi += manhattan_distance(x.second->get_v(), y.second->get_v());
      } //telos loop twn y
      mesi_apostasi /= (clust_points->size()-1);
      if (mesi_apostasi < min1)//min (mean(dist)) => new kentro
      {
        min1 = mesi_apostasi;
        min_dist_point.set_id(x.second->get_id());
        min_dist_point.set_v(x.second->get_v());
      }
    } //telos loop twn x
    //push back to neo kentro
    nea_kentra.push_back(min_dist_point);
  } //telos loop twn clusters
  format_clusters(&nea_kentra, clusts);
} //telos sunartisis


template <typename T>
void update_mean(std::vector<cluster<T>> * clusts, int diastaseis){

  std::vector<my_vector<T>> nea_kentra;
  for(unsigned int i=0; i< (*clusts).size(); i++){ //gia ka8e cluster
    T mean_coords[diastaseis] = {0.0};
    unsigned int n = (*clusts)[i].get_set_of_points()->size();
    std::unordered_map<std::string, my_vector<T> * > * clust_points = (*clusts)[i].get_set_of_points();
    for(auto x:*clust_points ){
      //std::cout << "henlo  " << x.second->get_v().size();
      for(unsigned int j=0; j< x.second->get_v().size(); j++ ){
        mean_coords[j] += x.second->get_v()[j] ;
      }//telos for gia tis suntetagmenes enos shmeiou tou cluster
      //n++;
    } //telos for gia kathe shmeiou tou cluster
    //std::cout << "eimai to n " << n << "\n";
    for(unsigned int z=0; z< diastaseis; z++)
      mean_coords[z] /= n;
    std::vector<T> mean_coordz(mean_coords, mean_coords+diastaseis); //metatroph se vector
    my_vector<T> mean_center;
    mean_center.set_id("\t");
    mean_center.set_v(mean_coordz);
    nea_kentra.push_back(mean_center);
  }//telos for gia ka8e cluster

  format_clusters(&nea_kentra, clusts);
  //return nea_kentra;
}//telos sunarthshs



////////////////CURVES/////////////////
template <typename T>
void update_pam_curve(std::vector<curve_cluster<T>> * clusts){
  std::vector<curve<T>> nea_kentra;
  for(unsigned int i=0; i< (*clusts).size(); i++){ //gia ka8e cluster
    double min1 = std::numeric_limits<double>::max(); //min pairnei timh apeiro
    curve<T> min_dist_point;
    std::unordered_map<std::string, curve<T> * > * clust_points = (*clusts)[i].get_set_of_curves(); //ta stoixeia tou cluster
    for(auto x:*clust_points) { //gia kathe stoixeio tou cluster
      double mesi_apostasi = 0.0;
      //sugkrine gia kathe shmeio != kentro
      if (x.first == (*clusts)[i].get_center_id())
        continue; //den psaxnw to idi kentro
      //else
      for(auto y:*clust_points) { //gia kathe stoixeio tou cluster
        //skip myself
        if (y.first == x.first)
          continue;
        //bres mean tou dist apo ta alla shmeia
        mesi_apostasi += dtw(x.second, y.second);
      } //telos loop twn y
      mesi_apostasi /= (clust_points->size()-1);
      if (mesi_apostasi < min1)//min (mean(dist)) => new kentro
      {
        min1 = mesi_apostasi;
        min_dist_point.set_id(x.second->get_id());
        min_dist_point.set_points(x.second->get_points());
      }
    } //telos loop twn x
    //push back to neo kentro
    nea_kentra.push_back(min_dist_point);
  } //telos loop twn clusters
  format_curve_clusters(&nea_kentra, clusts);
} //telos sunartisis
