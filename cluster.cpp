#include <iostream>
#include <string.h>
#include <fstream>
#include <typeinfo>
#include <stdlib.h> //rand?
#include <algorithm> // std::count
#include <cstring>
#include <limits>
#include <sstream>
#include <unordered_map>
#include <iomanip>
#include <set>
#include <cmath>
#include "my_vector.h"
#include "curve.h"
#include "cluster_object.h"
#include "utils.h"
#include "ht.h" //gia range search lsh
#include <chrono> // time measurements
#include <time.h>
#include "update.hpp"
#include "inits.hpp"
#include "assign.hpp"

std::string repeat_answer = "n";



bool is_ok_to_stop(std::vector<double> * objectives, int ith_rep){

  if(ith_rep >= 15)
    return true; //tha stamataei otan exei ftasei stis 25 epanalhpseis

  if(ith_rep < 2)
    return false;

  /*if(objectives->size() == 1)
    return false;

  if(objectives->size() == 2){
    if((*objectives)[0] < (*objectives)[1] )
      return false;
  }

  if(objectives->size() > 2){
    for(int i = objectives->size() -2; i>=0; i--){
      if( (*objectives)[i] == (*objectives)[objectives->size()-1] ){ //bre8hke ksana auth h timh (h teleutaia)
        if( (*objectives)[i+1] >= (*objectives)[objectives->size()-1] )
          return true; //tha ksanasunanthsoume th mikroterh timh kai tha stamathsoume tote
        else
          return false;
      }
    }
  }*/
  if((*objectives)[objectives->size()-2] - (*objectives)[objectives->size()-1] < 0)
    return false;

  if((*objectives)[objectives->size()-2] - (*objectives)[objectives->size()-1] < 0.01 * (*objectives)[objectives->size()-2]){
    return true;
  }

  return false;
}

int main(int argc, char *argv[])
{

  //test chamber

  //
  std::setprecision(14);
  bool complete = false;
  bool iset = false; ////an oxi orisma grammis entolos, 8a parw ta files apo path pou grafei o user
  bool cset = false;
  bool oset = false;
  char dataset_path[256];
  char config_path[256];
  char output_path[256];
  for (int i = 0; i < argc; i++)
  {
    if (strcmp("-i", argv[i]) == 0)
    {
      strcpy(dataset_path, argv[i + 1]);
      iset = true; /*std::cout <<"toxw\n";*/
    }              //pairnw to swsto onoma arxeiou apo to command line

    if (strcmp("-c", argv[i]) == 0)
    {
      strcpy(config_path, argv[i + 1]);
      cset = true; /*std::cout <<"toxw\n";*/
    }              //dinw timh sto K bash command line

    if (strcmp("-o", argv[i]) == 0)
    {
      strcpy(output_path, argv[i + 1]);
      oset = true; /*std::cout <<"toxw\n";*/
    }

    if (strcmp("-complete", argv[i]) == 0) //plithos hi functions gia dimiourgia twn g
    {
      complete = true;
    }
  }

  int number_of_grids = 2;
  int number_of_clusters = -1; //KALO EINAI NA YPARXEI STO ARXEIO
  int number_of_vector_hash_tables = 3;
  int number_of_vector_hash_functions = -1; //KALO EINAI NA YPARXEI STO ARXEIO
  if(cset){
    std::ifstream confile(config_path); //dataset: me tabs anamesa, ka8e grammi: id1    x11     x12     x13...
    std::string line;
    while (std::getline(confile, line))
    { //read files
      std::vector<std::string> tokens;
      std::stringstream check1(line);
      std::string intermediate;
      while (getline(check1, intermediate, ':'))
      {
          tokens.push_back(intermediate); //tokens has token[0] = id, token[1] = megethos, token[2] = all points (space separated now)
      }
      if(tokens[0] == "number_of_clusters")
        number_of_clusters = atoi(tokens[1].c_str());
      if(tokens[0] == "number_of_grids")
        number_of_grids = atoi(tokens[1].c_str());
      if(tokens[0] == "number_of_vector_hash_tables")
        number_of_vector_hash_tables =  atoi(tokens[1].c_str());
      if(tokens[0] == "number_of_vector_hash_functions")
        number_of_vector_hash_functions =  atoi(tokens[1].c_str());
    };
    confile.close();
    //ask values of essential vars (not default)
    while(number_of_vector_hash_functions < 0){
      std::cout << "number of hash functions not defined. Provide it:\n";
      std::cin >> number_of_vector_hash_functions;
    }
    while(number_of_clusters < 0){
      std::cout << "number of clusters not defined. Provide it:\n";
      std::cin >> number_of_clusters;
    }
  }
  else{ //keep default values and ask the 2 essentials (not default)
    while(number_of_vector_hash_functions < 0){
      std::cout << "number of hash functions not defined. Provide it:\n";
      std::cin >> number_of_vector_hash_functions;
    }
    while(number_of_clusters < 0){
      std::cout << "number of clusters not defined. Provide it:\n";
      std::cin >> number_of_clusters;
    }
  }

  //EAN DEN ORISTHKE APO GRAMMH ENTOLWN, DWSE MONOPATI DATASET:
  if (iset == false)
  {
    std::cout << "Define dataset path:\n";
    std::string inp1;
    std::cin >> dataset_path;
  }

  std::string what_is_the_input;
  int n = 0;                          //plithos twn vectors tou input file
  int diastaseis_vecs; //self explanatory
  std::ifstream infile(dataset_path); //dataset: me tabs anamesa, ka8e grammi: id1    x11     x12     x13...
  std::string line;
  bool datatype_set = false;
  //std::vector<my_vector<double>> vectors_array; //pinakas gia vectors
  std::unordered_map<std::string, my_vector<double> > vectors_array; //key == id tou vector, mapped value == my_vectors. Ta ids einai o,ti nai nai g auto....
  std::unordered_map<std::string, curve<double> > curves_array; //pinakas gia kampyles
  while (std::getline(infile, line))
  { //read files
    if((datatype_set == false)&&(line.find("vectors") != std::string::npos)){
      what_is_the_input = "vectors";
      datatype_set = true;
      continue;
    }

    //std::cout << line << "\n";
    if((datatype_set == false)&&(line.find("curves") != std::string::npos)){
      what_is_the_input = "curves";
      datatype_set = true;
      continue;
    }

    if(what_is_the_input == "vectors"){ // exoume na kanoyme me vectors
      my_vector<double> one_v_atime(line);
      //std::cout << one_v_atime.get_id()  <<"\n" ;
      diastaseis_vecs = one_v_atime.get_v().size();
      vectors_array[one_v_atime.get_id()] = one_v_atime;

      //vectors_array.push_back(one_v_atime);
      //std::cout << one_v_atime.get_id() << "\n";
      n++;
    }

    if(what_is_the_input == "curves"){ // exoume na kanoyme me vectors
      curve<double> one_v_atime(line);
      //std::cout << one_v_atime.get_id()  <<"\n" ;
      curves_array[one_v_atime.get_id()] = one_v_atime;
      n++;
    }

  };
  infile.close();
  //KSEKINAME ANALOGWS TO INPUT TYPE
  if(what_is_the_input == "vectors"){ //EXOUME NA KANOYME ME VECTORS

    std::vector<cluster<double>> clusters; //ta arxika mas clusters
    //INIT-1, ASS-1, UPD-1
    auto start = std::chrono::high_resolution_clock::now();
    initialise_centers<double>(number_of_clusters, &vectors_array, &clusters);
    double objective1 = lloyd_ass(&clusters, &vectors_array);
    std::vector<double> objectives;
    objectives.clear();
    objectives.push_back(objective1);
    std::cout << objective1 << "\n";
    int jot = 0;
    do{
      jot++;
      update_pam(&clusters);
      objective1 = lloyd_ass(&clusters, &vectors_array);
      objectives.push_back(objective1);
      std::cout << objective1 << "\n";
    }while(!(is_ok_to_stop(&objectives, jot)));
    auto end = std::chrono::high_resolution_clock::now() - start;
    long long microseconds = std::chrono::duration_cast<std::chrono::microseconds>(end).count();
    std::vector<double> sis = Silhouette(&clusters); //epistrefei ena vector me tis times twn SI gia kathe cluster
    //ektypwsh
    if (oset == false)
    {
      std::cout << "Define output file path:\n";
      std::string inp1;
      std::cin >> output_path;
    }
    std::ofstream outfile;
    outfile.open(output_path);
    outfile << "Algorithm: I1A1U1\n";
    for (unsigned int z = 0; z < clusters.size(); z++)
    {
      outfile << "CLUSTER-" << z << "{size: " << clusters[z].get_set_of_points()->size() << ", centroid: " << clusters[z].get_center_id() << "}"<<'\n';
    }
    outfile << "Clustering time: " << end.count() << '\n';
    outfile << "Silhouette: [";
    for (unsigned int z = 0; z < sis.size(); z++)
    {
      outfile << sis[z] << ",";
    }
    double sis_all = Silhouette_oliko<double>(sis);
    outfile << sis_all << "]\n";

    //optional
    if (complete == true)
    {
      for (unsigned int z = 0; z < clusters.size(); z++)
      {
        outfile << "CLUSTER-" << z << "{";
        for (auto x : *(clusters[z].get_set_of_points()))
        {
          outfile <<  x.second->get_id() << ", ";
        }
        outfile << "}\n";
      }
    }
    outfile << "\n";

    //INIT-1, ASS-1, UPD-2
    start = std::chrono::high_resolution_clock::now();
    initialise_centers<double>(number_of_clusters, &vectors_array, &clusters);
    objective1 = lloyd_ass(&clusters, &vectors_array);
    objectives.clear();
    objectives.push_back(objective1);
    std::cout << objective1 << "\n";
    jot = 0;
    do{
      jot++;
      update_mean(&clusters, diastaseis_vecs);
      objective1 = lloyd_ass(&clusters, &vectors_array);
      objectives.push_back(objective1);
      std::cout << objective1 << "\n";
    }while(!(is_ok_to_stop(&objectives, jot)));
    end = std::chrono::high_resolution_clock::now() - start;
    microseconds = std::chrono::duration_cast<std::chrono::microseconds>(end).count();
    sis.clear();
    sis = Silhouette(&clusters); //epistrefei ena vector me tis times twn SI gia kathe cluster
    //ektypwsh
    outfile << "Algorithm: I1A1U2\n";
    for (unsigned int z = 0; z < clusters.size(); z++)
    {
      outfile << "CLUSTER-" << z << "{size: " << clusters[z].get_set_of_points()->size() << ", centroid: ";
      for (unsigned int m = 0; m < clusters[z].get_center_coords().size(); m++)
      {
        outfile << clusters[z].get_center_coords()[m] << ", ";
      }
      outfile << "\n";
    }
    outfile << "Clustering time: " << end.count() << '\n';
    outfile << "Silhouette: [";
    for (unsigned int z = 0; z < sis.size(); z++)
    {
      outfile << sis[z] << ",";
    }
    sis_all = Silhouette_oliko<double>(sis);
    outfile << sis_all << "]\n";

    //optional
    if (complete == true)
    {
      for (unsigned int z = 0; z < clusters.size(); z++)
      {
        outfile << "CLUSTER-" << z << "{";
        for (auto x : *(clusters[z].get_set_of_points()))
        {
          outfile <<  x.second->get_id() << ", ";
        }
        outfile << "}\n";
      }
    }
    outfile << "\n";


    //INIT-1, ASS-2, UPD-1
    start = std::chrono::high_resolution_clock::now();
    initialise_centers<double>(number_of_clusters, &vectors_array, &clusters);
    double w_lsh = - 1.0;
    objective1 = LSH_range_ass(&clusters, &vectors_array, diastaseis_vecs, number_of_vector_hash_tables, number_of_vector_hash_functions, &w_lsh);
    objectives.clear();
    objectives.push_back(objective1);
    std::cout << objective1 << "\n";
    jot = 0;
    do{
      jot++;
      update_pam(&clusters);
      objective1 = LSH_range_ass(&clusters, &vectors_array, diastaseis_vecs, number_of_vector_hash_tables, number_of_vector_hash_functions, &w_lsh);
      objectives.push_back(objective1);
      std::cout << objective1 << "\n";
    }while(!(is_ok_to_stop(&objectives, jot)));
    end = std::chrono::high_resolution_clock::now() - start;
    microseconds = std::chrono::duration_cast<std::chrono::microseconds>(end).count();
    sis.clear();
    sis = Silhouette(&clusters); //epistrefei ena vector me tis times twn SI gia kathe cluster
    //ektypwsh
    outfile << "Algorithm: I1A2U1\n";
    for (unsigned int z = 0; z < clusters.size(); z++)
    {
      outfile << "CLUSTER-" << z << "{size: " << clusters[z].get_set_of_points()->size() << ", centroid: " << clusters[z].get_center_id() << "}"<<'\n';
    }
    outfile << "Clustering time: " << end.count() << '\n';
    outfile << "Silhouette: [";
    for (unsigned int z = 0; z < sis.size(); z++)
    {
      outfile << sis[z] << ",";
    }
    sis_all = Silhouette_oliko<double>(sis);
    outfile << sis_all << "]\n";

    //optional
    if (complete == true)
    {
      for (unsigned int z = 0; z < clusters.size(); z++)
      {
        outfile << "CLUSTER-" << z << "{";
        for (auto x : *(clusters[z].get_set_of_points()))
        {
          outfile <<  x.second->get_id() << ", ";
        }
        outfile << "}\n";
      }
    }
    outfile << "\n";


    //INIT-1, ASS-2, UPD-2
    start = std::chrono::high_resolution_clock::now();
    initialise_centers<double>(number_of_clusters, &vectors_array, &clusters);
    //w_lsh = - 1.0;
    objective1 = LSH_range_ass(&clusters, &vectors_array, diastaseis_vecs, number_of_vector_hash_tables, number_of_vector_hash_functions, &w_lsh);
    objectives.clear();
    objectives.push_back(objective1);
    std::cout << objective1 << "\n";
    jot = 0;
    do{
      jot++;
      update_mean(&clusters, diastaseis_vecs);
      objective1 = LSH_range_ass(&clusters, &vectors_array, diastaseis_vecs, number_of_vector_hash_tables, number_of_vector_hash_functions, &w_lsh);
      objectives.push_back(objective1);
      std::cout << objective1 << "\n";
    }while(!(is_ok_to_stop(&objectives, jot)));
    end = std::chrono::high_resolution_clock::now() - start;
    microseconds = std::chrono::duration_cast<std::chrono::microseconds>(end).count();
    sis.clear();
    sis = Silhouette(&clusters); //epistrefei ena vector me tis times twn SI gia kathe cluster
    //ektypwsh
    outfile << "Algorithm: I1A2U2\n";
    for (unsigned int z = 0; z < clusters.size(); z++)
    {
      outfile << "CLUSTER-" << z << "{size: " << clusters[z].get_set_of_points()->size() << ", centroid: ";
      for (unsigned int m = 0; m < clusters[z].get_center_coords().size(); m++)
      {
        outfile << clusters[z].get_center_coords()[m] << ", ";
      }
      outfile << "\n";
    }
    outfile << "Clustering time: " << end.count() << '\n';
    outfile << "Silhouette: [";
    for (unsigned int z = 0; z < sis.size(); z++)
    {
      outfile << sis[z] << ",";
    }
    sis_all = Silhouette_oliko<double>(sis);
    outfile << sis_all << "]\n";

    //optional
    if (complete == true)
    {
      for (unsigned int z = 0; z < clusters.size(); z++)
      {
        outfile << "CLUSTER-" << z << "{";
        for (auto x : *(clusters[z].get_set_of_points()))
        {
          outfile <<  x.second->get_id() << ", ";
        }
        outfile << "}\n";
      }
    }
    outfile << "\n";


    //INIT-2, ASS-1, UPD-1
    start = std::chrono::high_resolution_clock::now();
    initialise_centers_plus<double>(number_of_clusters, &vectors_array, &clusters);
    objective1 = lloyd_ass(&clusters, &vectors_array);
    objectives.clear();
    objectives.push_back(objective1);
    std::cout << objective1 << "\n";
    jot = 0;
    do{
      jot++;
      update_pam(&clusters);
      objective1 = lloyd_ass(&clusters, &vectors_array);
      objectives.push_back(objective1);
      std::cout << objective1 << "\n";
    }while(!(is_ok_to_stop(&objectives, jot)));
    end = std::chrono::high_resolution_clock::now() - start;
    microseconds = std::chrono::duration_cast<std::chrono::microseconds>(end).count();
    sis.clear();
    sis = Silhouette(&clusters); //epistrefei ena vector me tis times twn SI gia kathe cluster
    //ektypwsh
    outfile << "Algorithm: I2A1U1\n";
    for (unsigned int z = 0; z < clusters.size(); z++)
    {
      outfile << "CLUSTER-" << z << "{size: " << clusters[z].get_set_of_points()->size() << ", centroid: " << clusters[z].get_center_id()<<"}" <<'\n';
    }
    outfile << "Clustering time: " << end.count() << '\n';
    outfile << "Silhouette: [";
    for (unsigned int z = 0; z < sis.size(); z++)
    {
      outfile << sis[z] << ",";
    }
    sis_all = Silhouette_oliko<double>(sis);
    outfile << sis_all << "]\n";

    //optional
    if (complete == true)
    {
      for (unsigned int z = 0; z < clusters.size(); z++)
      {
        outfile << "CLUSTER-" << z << "{";
        for (auto x : *(clusters[z].get_set_of_points()))
        {
          outfile <<  x.second->get_id() << ", ";
        }
        outfile << "}\n";
      }
    }
    outfile << "\n";


    //INIT-2, ASS-1, UPD-2
    start = std::chrono::high_resolution_clock::now();
    initialise_centers_plus<double>(number_of_clusters, &vectors_array, &clusters);
    objective1 = lloyd_ass(&clusters, &vectors_array);
    objectives.clear();
    objectives.push_back(objective1);
    std::cout << objective1 << "\n";
    jot = 0;
    do{
      jot++;
      update_mean(&clusters, diastaseis_vecs);
      objective1 = lloyd_ass(&clusters, &vectors_array);
      objectives.push_back(objective1);
      std::cout << objective1 << "\n";
    }while(!(is_ok_to_stop(&objectives, jot)));
    end = std::chrono::high_resolution_clock::now() - start;
    microseconds = std::chrono::duration_cast<std::chrono::microseconds>(end).count();
    sis.clear();
    sis = Silhouette(&clusters); //epistrefei ena vector me tis times twn SI gia kathe cluster
    //ektypwsh
    outfile << "Algorithm: I2A1U2\n";
    for (unsigned int z = 0; z < clusters.size(); z++)
    {
      outfile << "CLUSTER-" << z << "{size: " << clusters[z].get_set_of_points()->size() << ", centroid: ";
      for (unsigned int m = 0; m < clusters[z].get_center_coords().size(); m++)
      {
        outfile << clusters[z].get_center_coords()[m] << ", ";
      }
      outfile << "\n";
    }
    outfile << "Clustering time: " << end.count() << '\n';
    outfile << "Silhouette: [";
    for (unsigned int z = 0; z < sis.size(); z++)
    {
      outfile << sis[z] << ",";
    }
    sis_all = Silhouette_oliko<double>(sis);
    outfile << sis_all << "]\n";

    //optional
    if (complete == true)
    {
      for (unsigned int z = 0; z < clusters.size(); z++)
      {
        outfile << "CLUSTER-" << z << "{";
        for (auto x : *(clusters[z].get_set_of_points()))
        {
          outfile <<  x.second->get_id() << ", ";
        }
        outfile << "}\n";
      }
    }
    outfile << "\n";


    //INIT-2, ASS-2, UPD-1
    start = std::chrono::high_resolution_clock::now();
    initialise_centers_plus<double>(number_of_clusters, &vectors_array, &clusters);
    //w_lsh = - 1.0;
    objective1 = LSH_range_ass(&clusters, &vectors_array, diastaseis_vecs, number_of_vector_hash_tables, number_of_vector_hash_functions, &w_lsh);
    objectives.clear();
    objectives.push_back(objective1);
    std::cout << objective1 << "\n";
    jot = 0;
    do{
      jot++;
      update_pam(&clusters);
      objective1 = LSH_range_ass(&clusters, &vectors_array, diastaseis_vecs, number_of_vector_hash_tables, number_of_vector_hash_functions, &w_lsh);
      objectives.push_back(objective1);
      std::cout << objective1 << "\n";
    }while(!(is_ok_to_stop(&objectives, jot)));
    end = std::chrono::high_resolution_clock::now() - start;
    microseconds = std::chrono::duration_cast<std::chrono::microseconds>(end).count();
    sis.clear();
    sis = Silhouette(&clusters); //epistrefei ena vector me tis times twn SI gia kathe cluster
    //ektypwsh
    outfile << "Algorithm: I2A2U1\n";
    for (unsigned int z = 0; z < clusters.size(); z++)
    {
      outfile << "CLUSTER-" << z << "{size: " << clusters[z].get_set_of_points()->size() << ", centroid: " << clusters[z].get_center_id() << "}" <<'\n';
    }
    outfile << "Clustering time: " << end.count() << '\n';
    outfile << "Silhouette: [";
    for (unsigned int z = 0; z < sis.size(); z++)
    {
      outfile << sis[z] << ",";
    }
    sis_all = Silhouette_oliko<double>(sis);
    outfile << sis_all << "]\n";

    //optional
    if (complete == true)
    {
      for (unsigned int z = 0; z < clusters.size(); z++)
      {
        outfile << "CLUSTER-" << z << "{";
        for (auto x : *(clusters[z].get_set_of_points()))
        {
          outfile <<  x.second->get_id() << ", ";
        }
        outfile << "}\n";
      }
    }
    outfile << "\n";


    //INIT-2, ASS-2, UPD-2
    start = std::chrono::high_resolution_clock::now();
    initialise_centers_plus<double>(number_of_clusters, &vectors_array, &clusters);
    //w_lsh = - 1.0;
    objective1 = LSH_range_ass(&clusters, &vectors_array, diastaseis_vecs, number_of_vector_hash_tables, number_of_vector_hash_functions, &w_lsh);
    objectives.clear();
    objectives.push_back(objective1);
    std::cout << objective1 << "\n";
    jot = 0;
    do{
      jot++;
      update_mean(&clusters, diastaseis_vecs);
      objective1 = LSH_range_ass(&clusters, &vectors_array, diastaseis_vecs, number_of_vector_hash_tables, number_of_vector_hash_functions, &w_lsh);
      objectives.push_back(objective1);
      std::cout << objective1 << "\n";
    }while(!(is_ok_to_stop(&objectives, jot)));
    end = std::chrono::high_resolution_clock::now() - start;
    microseconds = std::chrono::duration_cast<std::chrono::microseconds>(end).count();
    sis.clear();
    sis = Silhouette(&clusters); //epistrefei ena vector me tis times twn SI gia kathe cluster
    //ektypwsh
    outfile << "Algorithm: I2A2U2\n";
    for (unsigned int z = 0; z < clusters.size(); z++)
    {
      outfile << "CLUSTER-" << z << "{size: " << clusters[z].get_set_of_points()->size() << ", centroid: ";
      for (unsigned int m = 0; m < clusters[z].get_center_coords().size(); m++)
      {
        outfile << clusters[z].get_center_coords()[m]<<", ";
      }
      outfile << "\n";
    }
    outfile << "Clustering time: " << end.count() << '\n';
    outfile << "Silhouette: [";
    for (unsigned int z = 0; z < sis.size(); z++)
    {
      outfile << sis[z] << ",";
    }
    sis_all = Silhouette_oliko<double>(sis);
    outfile << sis_all << "]\n";

    //optional
    if (complete == true)
    {
      for (unsigned int z = 0; z < clusters.size(); z++)
      {
        outfile << "CLUSTER-" << z << "{";
        for (auto x : *(clusters[z].get_set_of_points()))
        {
          outfile <<  x.second->get_id() << ", ";
        }
        outfile << "}\n";
      }
    }
    outfile << "\n";

    /*for(unsigned int i = 0; i < clusters.size(); i++){
      clusters[i].print_cluster();
    }*/
    outfile.close();
  }
  else if(what_is_the_input == "curves"){ //an mas do8oun kampyles

    std::vector<curve_cluster<double>> clusters; //ta arxika mas clusters
    std::pair<double, double> twoelems = calculate_delta(&curves_array);
    double delta = twoelems.second;
    double max_coord_lsh = twoelems.first;
    double w_lsh = - 1.0;
    std::vector<double> objectives;

    //INIT-1, ASS-1, UPD-1
    auto start = std::chrono::high_resolution_clock::now();
    initialise_centers_curve(number_of_clusters, &curves_array, &clusters);
    double objective1 = lloyd_ass_curve(&clusters, &curves_array);
    objectives.clear();
    objectives.push_back(objective1);
    std::cout << objective1 << "\n";
    int jot = 0;
    do{
      jot++;
      update_pam_curve(&clusters);
      objective1 =  lloyd_ass_curve(&clusters, &curves_array);
      //objective1 = LSH_range_ass_curve(&clusters, &curves_array, number_of_grids, number_of_vector_hash_functions, delta, max_coord_lsh, &w_lsh);
      objectives.push_back(objective1);
      std::cout << objective1 << "\n";
    }while(!(is_ok_to_stop(&objectives, jot)));
    auto end = std::chrono::high_resolution_clock::now() - start;
    long long microseconds = std::chrono::duration_cast<std::chrono::microseconds>(end).count();
    std::vector<double> sis = Silhouette_curve(&clusters); //epistrefei ena vector me tis times twn SI gia kathe cluster
    //ektypwsh
    if (oset == false)
    {
      std::cout << "Define output file path:\n";
      std::string inp1;
      std::cin >> output_path;
    }
    std::ofstream outfile;
    outfile.open(output_path);
    outfile << "Algorithm: I1A1U1\n";
    //ektypwsh
    for (unsigned int z = 0; z < clusters.size(); z++)
    {
      outfile << "CLUSTER-" << z << "{size: " << clusters[z].get_set_of_curves()->size() << ", centroid: " << clusters[z].get_center_id() << "}" <<'\n';
    }
    outfile << "Clustering time: " << end.count() << '\n';
    outfile << "Silhouette: [";
    for (unsigned int z = 0; z < sis.size(); z++)
    {
      outfile << sis[z] << ",";
    }
    double sis_all = Silhouette_oliko<double>(sis);
    outfile << sis_all << "]\n";

    //optional
    if (complete == true)
    {
      for (unsigned int z = 0; z < clusters.size(); z++)
      {
        outfile << "CLUSTER-" << z << "{";
        for (auto x : *(clusters[z].get_set_of_curves()))
        {
          outfile <<  x.second->get_id() << ", ";
        }
        outfile << "}\n";
      }
    }
    outfile << "\n";


    //INIT-1, ASS-2, UPD-1
    start = std::chrono::high_resolution_clock::now();
    initialise_centers_curve(number_of_clusters, &curves_array, &clusters);
    w_lsh = - 1.0;
    objective1 = LSH_range_ass_curve(&clusters, &curves_array, number_of_grids, number_of_vector_hash_functions, delta, max_coord_lsh, &w_lsh);
    objectives.clear();
    objectives.push_back(objective1);
    std::cout << objective1 << "\n";
    jot = 0;
    do{
      jot++;
      update_pam_curve(&clusters);
      objective1 = LSH_range_ass_curve(&clusters, &curves_array, number_of_grids, number_of_vector_hash_functions, delta, max_coord_lsh, &w_lsh);
      objectives.push_back(objective1);
      std::cout << objective1 << "\n";
    }while(!(is_ok_to_stop(&objectives, jot)));
    end = std::chrono::high_resolution_clock::now() - start;
    microseconds = std::chrono::duration_cast<std::chrono::microseconds>(end).count();
    sis.clear();
    sis = Silhouette_curve(&clusters); //epistrefei ena vector me tis times twn SI gia kathe cluster
    outfile << "Algorithm: I1A2U1\n";
    for (unsigned int z = 0; z < clusters.size(); z++)
    {
      outfile << "CLUSTER-" << z << "{size: " << clusters[z].get_set_of_curves()->size() << ", centroid: " << clusters[z].get_center_id() << "}" <<'\n';
    }
    outfile << "Clustering time: " << end.count() << '\n';
    outfile << "Silhouette: [";
    for (unsigned int z = 0; z < sis.size(); z++)
    {
      outfile << sis[z] << ",";
    }
    sis_all = Silhouette_oliko<double>(sis);
    outfile << sis_all << "]\n";
    //optional
    if (complete == true)
    {
      for (unsigned int z = 0; z < clusters.size(); z++)
      {
        outfile << "CLUSTER-" << z << "{";
        for (auto x : *(clusters[z].get_set_of_curves()))
        {
          outfile <<  x.second->get_id() << ", ";
        }
        outfile << "}\n";
      }
    }
    outfile << "\n";


    //INIT-2, ASS-1, UPD-1
    start = std::chrono::high_resolution_clock::now();
    initialise_centers_plus_curve(number_of_clusters, &curves_array, &clusters);
    objective1 = lloyd_ass_curve(&clusters, &curves_array);
    objectives.clear();
    objectives.push_back(objective1);
    std::cout << objective1 << "\n";
    jot = 0;
    do{
      jot++;
      update_pam_curve(&clusters);
      objective1 = lloyd_ass_curve(&clusters, &curves_array);
      objectives.push_back(objective1);
      std::cout << objective1 << "\n";
    }while(!(is_ok_to_stop(&objectives, jot)));
    end = std::chrono::high_resolution_clock::now() - start;
    microseconds = std::chrono::duration_cast<std::chrono::microseconds>(end).count();
    sis.clear();
    sis = Silhouette_curve(&clusters); //epistrefei ena vector me tis times twn SI gia kathe cluster
    outfile << "Algorithm: I2A1U1\n";
    for (unsigned int z = 0; z < clusters.size(); z++)
    {
      outfile << "CLUSTER-" << z << "{size: " << clusters[z].get_set_of_curves()->size() << ", centroid: " << clusters[z].get_center_id() << "}" <<'\n';
    }
    outfile << "Clustering time: " << end.count() << '\n';
    outfile << "Silhouette: [";
    for (unsigned int z = 0; z < sis.size(); z++)
    {
      outfile << sis[z] << ",";
    }
    sis_all = Silhouette_oliko<double>(sis);
    outfile << sis_all << "]\n";
    //optional
    if (complete == true)
    {
      for (unsigned int z = 0; z < clusters.size(); z++)
      {
        outfile << "CLUSTER-" << z << "{";
        for (auto x : *(clusters[z].get_set_of_curves()))
        {
          outfile <<  x.second->get_id() << ", ";
        }
        outfile << "}\n";
      }
    }
    outfile << "\n";


    //INIT-2, ASS-2, UPD-1
    start = std::chrono::high_resolution_clock::now();
    initialise_centers_plus_curve(number_of_clusters, &curves_array, &clusters);
    //w_lsh = - 1.0;
    objective1 = LSH_range_ass_curve(&clusters, &curves_array, number_of_grids, number_of_vector_hash_functions, delta, max_coord_lsh, &w_lsh);
    objectives.clear();
    objectives.push_back(objective1);
    std::cout << objective1 << "\n";
    jot = 0;
    do{
      jot++;
      update_pam_curve(&clusters);
      objective1 = LSH_range_ass_curve(&clusters, &curves_array, number_of_grids, number_of_vector_hash_functions, delta, max_coord_lsh, &w_lsh);
      objectives.push_back(objective1);
      std::cout << objective1 << "\n";
    }while(!(is_ok_to_stop(&objectives, jot)));
    end = std::chrono::high_resolution_clock::now() - start;
    microseconds = std::chrono::duration_cast<std::chrono::microseconds>(end).count();
    sis.clear();
    sis = Silhouette_curve(&clusters); //epistrefei ena vector me tis times twn SI gia kathe cluster
    outfile << "Algorithm: I2A2U1\n";
    for (unsigned int z = 0; z < clusters.size(); z++)
    {
      outfile << "CLUSTER-" << z << "{size: " << clusters[z].get_set_of_curves()->size() << ", centroid: " << clusters[z].get_center_id() << "}" <<'\n';
    }
    outfile << "Clustering time: " << end.count() << '\n';
    outfile << "Silhouette: [";
    for (unsigned int z = 0; z < sis.size(); z++)
    {
      outfile << sis[z] << ",";
    }
    sis_all = Silhouette_oliko<double>(sis);
    outfile << sis_all << "]\n";
    //optional
    if (complete == true)
    {
      for (unsigned int z = 0; z < clusters.size(); z++)
      {
        outfile << "CLUSTER-" << z << "{";
        for (auto x : *(clusters[z].get_set_of_curves()))
        {
          outfile <<  x.second->get_id() << ", ";
        }
        outfile << "}\n";
      }
    }
    outfile << "\n";

    /*for(unsigned int i = 0; i < clusters.size(); i++){
      clusters[i].print_cluster();
    }*/
  }
  else{
    std::cout << "Den orises ti typou dedomena exoyme sthn prwth grammh opws eipe h ekfnwhsh, Enjoy the exit :* xoxo\n";
    exit(-1);
  }

}
