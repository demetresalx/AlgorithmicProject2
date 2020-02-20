#include "utils.h"
#include <unordered_map>
#include "cluster_object.h"
#include "my_vector.h"
#include "NNpair.h"
#include "ht.h"
#include "curve_ht.h"
#include "grid.h"


//LLOYD - each input vector to it nearest!
template <typename T>
double lloyd_ass(std::vector<cluster<T>>* clusters, std::unordered_map<std::string, my_vector<T> > *vectors_array)
{
    double obj_fun = 0.0;
    for (auto x :(*vectors_array))
    {
        double min1 = std::numeric_limits<double>::max(); //apeiro
        int min_clust_index = -1; //to index tou cluster opou anoikei kathe x
        for (unsigned int j = 0; j < (*clusters).size(); j++)
        {
            double tmp = 0.0;
            tmp = manhattan_distance(x.second.get_v(), (*clusters)[j].get_center_coords()); //dist(x.second.get_v(), (*clusters)[j].get_center_coords())
            /*std::vector<T> weko = (*clusters)[j].get_center_coords();
            for(unsigned int haha=0; haha< weko.size(); haha++){
              std::cout << weko[haha] << " ";
            }
            std::cout << "\n";*/
            if (tmp < min1)
            {
                min1 = tmp;
                min_clust_index = j;
            }
            //std::cout << "Eimai " << x.first << " & sugkrinw " << (*clusters)[j].get_center_id() << " dist " << tmp << '\n';
        }
        (*clusters)[min_clust_index].incorporate_point( &((*vectors_array)[x.first]) );
        obj_fun += min1;
        //std::cout << "Eimai " << x.first << " & kentro " << (*clusters)[min_clust_index].get_center_id() << " dist " << min1 << '\n';
        //std::cout << '\n';
    }
    //std::cout << "\n\n\n";
    int biggest_index = 0;
    int maxim_numb = (*clusters)[0].get_set_of_points()->size();
    for (unsigned int j = 1; j < (*clusters).size(); j++){
      if((*clusters)[j].get_set_of_points()->size() > maxim_numb){
        maxim_numb = (*clusters)[j].get_set_of_points()->size();
        biggest_index = j;
      }
    }

    //prepei na sigoureutw oti den yparxei adeio cluster. An yparxei, pare to makrynotero shmeio tou pio megalou cluster. Einai apodekth methodos
    for (unsigned int j = 0; j < (*clusters).size(); j++){
      std::unordered_map<std::string, my_vector<T> * > * clust_points = (*clusters)[j].get_set_of_points();
      if((* clust_points).size() < 1){
        //std::cout << "yparxei adeio cluster!\n";
        std::unordered_map<std::string, my_vector<T> * > * big_points = (*clusters)[biggest_index].get_set_of_points();
        double max_distn = std::numeric_limits<double>::max();           //apeiro kai kala
        max_distn = -1 * max_distn;                               //arxika -apeiro
        std::string worst_id;
        for(auto xy: *big_points){
          double dis_dist = manhattan_distance(xy.second->get_v(), (*clusters)[biggest_index].get_center_coords());
          if(dis_dist > max_distn){
            max_distn = dis_dist;
            worst_id = xy.first;
          } //telos if apostashs
        } //telos for anazhthshs sto pio xontrouli cluster
        (*clusters)[biggest_index].discorporate_point(&((*vectors_array)[worst_id]));
        obj_fun -= max_distn;
        (*clusters)[j].incorporate_point(&((*vectors_array)[worst_id]));
        obj_fun += manhattan_distance((*vectors_array)[worst_id].get_v(), (*clusters)[j].get_center_coords());;
      }//telos if periptwshs adeiou cluster
    } //telos for gia anazhthsh adeiwn clusters
    return obj_fun;
}



//INVERSE ASSIGNMENT WITH RANGE SEARCH LSH BLA BLA
template <typename T> //briskei thn arxikh aktina
double initialize_radius(std::vector<cluster<T>>* clusters)
{
  double init_radius = std::numeric_limits<double>::max();
  for(unsigned int i=0; i< clusters->size(); i++){
    for(unsigned int j=0; j< clusters->size(); j++){

      if((*clusters)[i].get_center_id() == (*clusters)[j].get_center_id()) //einai to idio kentro
        continue;

      if( manhattan_distance((*clusters)[i].get_center_coords(), (*clusters)[j].get_center_coords()) < init_radius)
          init_radius = manhattan_distance((*clusters)[i].get_center_coords(), (*clusters)[j].get_center_coords()) ;

    }
  }
  return init_radius / 2;
}

//se curve edition
template <typename T> //briskei thn arxikh aktina
double initialize_radius_curve(std::vector<curve_cluster<T>>* clusters)
{
  double init_radius = std::numeric_limits<double>::max();
  for(unsigned int i=0; i< clusters->size(); i++){
    for(unsigned int j=0; j< clusters->size(); j++){

      if((*clusters)[i].get_center_id() == (*clusters)[j].get_center_id()) //einai to idio kentro
        continue;

      double apostasi= dtw((*clusters)[i].get_center_ptr(), (*clusters)[j].get_center_ptr());
      if( apostasi < init_radius)
          init_radius = apostasi ;

    }
  }
  return init_radius / 2;
}




template <typename T> //APO ERGASIA 1
double LSH_range_ass(std::vector<cluster<T>>* clusters, std::unordered_map<std::string, my_vector<T> > *vectors_array,int diastaseis ,int number_of_vector_hash_tables, int number_of_vector_hash_functions, double *w_done)
{
    double w = *w_done;
    if(*w_done < 0){ //den exei ypologistei to w, pame na to broume
      //auto start_of_w_calc = std::chrono::high_resolution_clock::now();
      std::vector<NNpair> input_actual_NNs; //pinakas apo zeugaria actual NNs me prwto stoixeio to p
      for (auto x :(*vectors_array))
      { //prepei na brw ta zeugaria ap to input gia ypologismo w

        std::string min_id1;
        double min1 = std::numeric_limits<double>::max(); //min pairnei timh apeiro
        for (auto y :(*vectors_array))
        {
          if (manhattan_distance(x.second.get_v(), y.second.get_v()) == 0) //einai to idio shmeio
            continue;

          if (manhattan_distance(x.second.get_v(), y.second.get_v()) < min1)
          {
            min1 = manhattan_distance(x.second.get_v(), y.second.get_v());
            min_id1 =  y.second.get_id();
          }
        }
        NNpair single_pair1(x.second.get_id(), min_id1);
        single_pair1.set_distance(min1);
        input_actual_NNs.push_back(single_pair1);
      }

      double tmp = 0.0;
      double mean_distance = 0;
      for (unsigned int i = 0; i < input_actual_NNs.size(); i++)
      {
        tmp += input_actual_NNs.at(i).get_distance();
      }
      //int diastaseis = vectors_array->begin().second.get_v().size(); //gia hash tables
      //std::cout << "diasteaseis " << diastaseis << "\n";

      mean_distance = tmp / input_actual_NNs.size(); //fp division
      //auto end_of_w_calc = std::chrono::high_resolution_clock::now() - start_of_w_calc;
      //long long microseconds_w = std::chrono::duration_cast<std::chrono::microseconds>(end_of_w_calc).count();
      //fprintf(stderr, "Time needed for w calculation is %lld microseconds.\n\n", microseconds_w);
      fprintf(stderr, "Value of w = %f\n", mean_distance);
      //also test gia w = 10 * mean_distance
      /*const*/ w = 4 * mean_distance; //to w pou vazw sta ai, STH XEIROTERH HARD CODED
      *w_done = w; //ypologisthke prwth kai teleutaia fora!
    }


    /////////////////////////////LSH TIME////////////////////////////////
    int Table_Size = (*vectors_array).size() / 8;
    std::vector<ht<T>> our_hash_tables;
    for (int i = 0; i < number_of_vector_hash_tables; i++)
    {
      ht<T> temp_hash_table(Table_Size, number_of_vector_hash_functions, diastaseis, w);
      our_hash_tables.push_back(temp_hash_table);
    }

    std::unordered_map<std::string, std::pair<int, double>> owned; //to flag poy lene oi diafaneies gia to an kaparw8hke ena shmeio kai apo poio index kai me poia aktina
    int n =0;
    //for (auto x :(*vectors_array)) //hasharw ta vectors k arxikopoiw flags
    for (auto x :(*vectors_array)) //hasharw ta vectors k arxikopoiw flags
    {
      for (int j = 0; j < number_of_vector_hash_tables; j++)
        {our_hash_tables[j].hash_vector(&((*vectors_array)[x.first]));}

      std::pair<int, double> index_and_radius;
      index_and_radius.first = -1;
      index_and_radius.second = 0.0;
      owned[x.first] = index_and_radius;
      n++;
    }

    int num_unassigned = n; //posa exoun meinei xwris anathesh
    double radius = initialize_radius(clusters);
    std::vector<std::string> this_center_neighbs;
    std::vector<std::string> this_HT_neighbs;
    bool repetition = false;
    int kill_countdown = 15; //an den exoun ginei nees anatheseis meta apo tosous diplasiasmous aktinas, stop
    int num_unassigned_prev = num_unassigned; //arithmos unassigned shmeiwn prin th loypa gia na sugkrinoyme proodo kathe fora kai na stamatame

    while((num_unassigned > n/10) && (kill_countdown >0) ){ //h anazhthsh range search lsh tha ginetai mexri to 90% twn shmeiwn ginei assign se kapoio kentro. Epeita klassikh methodos opws prota8hke

      num_unassigned_prev = num_unassigned;
      for(unsigned int i=0; i< clusters->size(); i++){ //gia kathe kentro twn clusters

        this_center_neighbs.clear();
        for (int j = 0; j < number_of_vector_hash_tables; j++){ //LSH se L hashtables
          this_HT_neighbs.clear();
          this_HT_neighbs = our_hash_tables[j].hash_query((*clusters)[i].get_center_ptr(), radius, repetition);
          //std::cout << this_HT_neighbs.size() << "-";
          this_center_neighbs.insert(this_center_neighbs.end(), this_HT_neighbs.begin(), this_HT_neighbs.end());
          //std::cout << this_center_neighbs.size() << " ";
        }
        //pros8hkh shmeiwn se cluster kai flag gia na mhn to paroyn kai ta ypoloipa clusters
        //std::cout << "eimai to cl " << (*clusters)[i].get_center_id() << "kai "<<this_center_neighbs.size() << "\n";
        for(unsigned int z=0; z< this_center_neighbs.size(); z++){
          //std::cout << owned[this_center_neighbs[z]].first ;
          //std::cout << this_center_neighbs[z] ;
          if( owned[this_center_neighbs[z]].first == -1 ){ //den exei kaparw8ei
              //std::cout << "kapakap ";
              (*clusters)[i].incorporate_point(&((*vectors_array)[this_center_neighbs[z]]));
              owned[this_center_neighbs[z]].first = i; //to kaparwse
              owned[this_center_neighbs[z]].second = radius; //to kaparwse entos aktinas toshs
              num_unassigned--;
          }
          else{ //sugkrinoume me auton poy to exei kaparwsei
            if(owned[this_center_neighbs[z]].second >= radius){ //an kapoios allos to exei kaparwsei me mikroterh aktina, apofeugoume th sugkrish kai proxwrame
              if(manhattan_distance((*vectors_array)[this_center_neighbs[z]].get_v() , (*clusters)[i].get_center_coords()) < manhattan_distance((*vectors_array)[this_center_neighbs[z]].get_v() , (*clusters)[owned[this_center_neighbs[z]].first].get_center_coords()) ){
                (*clusters)[owned[this_center_neighbs[z]].first].discorporate_point(&((*vectors_array)[this_center_neighbs[z]])); //to bgazei ap to palio
                (*clusters)[i].incorporate_point(&((*vectors_array)[this_center_neighbs[z]])); //to vazei sto neo
                owned[this_center_neighbs[z]].first = i; //to kaparwse
                owned[this_center_neighbs[z]].second = radius; //to kaparwse entos aktinas toshs
              } //telos if gia sugkrish apostasewn
            } //telos if poy afora an ena kaparwmeno shmeio exei kaparw8ei apo mikroterh aktina ara den exei nohma na koitaksoume pali
          } //telos else poy afora to an ena shmeio einai kaparwmeno h oxi
        } //telos gor gia auta poy brhke auto to cluster gia authn thn aktina
      } //telos for gia ta clusters

      if(num_unassigned_prev - num_unassigned <=0 ) //den kaname proodo, arxise antistrofh metrhsh
        kill_countdown--;
      else
        kill_countdown = 15; //eixame proodo, mhdenise thn antistrofh metrhsh

      if(kill_countdown<=0) //den yphrkse veltiwsh gia sunexomenes loypes, telos
        break;

      radius = radius*2; //diplasiazoume aktina kai sunexizoume
      repetition = true;
      //std::cout << num_unassigned << "\n";
    } //telos ths while poy diplasiazoume thn aktina


    //twra an kapoio exei meinei akaparwto, prepei na paei sto kontinotero tou kentro
    double objective_function = 0.0;
    for(auto x: owned){
      if(x.second.first == -1){ //akaparwto
        double min1 = std::numeric_limits<double>::max(); //apeiro
        int min_clust_index = -1; //to index tou cluster opou anoikei kathe x
        //(*vectors_array)[x.first].get_id()
        for (unsigned int j = 0; j < (*clusters).size(); j++)
        {
            double tmp = 0.0;
            tmp = manhattan_distance( (*vectors_array)[x.first].get_v(), (*clusters)[j].get_center_coords()); //dist(x.second.get_v(), (*clusters)[j].get_center_coords())
            if (tmp < min1)
            {
                min1 = tmp;
                min_clust_index = j;
            }
            //std::cout << "Eimai " << x.first << " & sugkrinw " << (*clusters)[j].get_center_id() << " dist " << tmp << '\n';
        }
        (*clusters)[min_clust_index].incorporate_point(&((*vectors_array)[x.first]));
        objective_function += min1;
      }
      else{ //kaparwmeno, ypologizw aplws apostash gia to objectve function
        objective_function += manhattan_distance( (*vectors_array)[x.first].get_v() , (*clusters)[x.second.first].get_center_coords() ); //h  apostash metaksu tou shmeiou kai tou kentro poy to kaparwse
      }
    }//telos for gia akaparwta

    int biggest_index = 0;
    int maxim_numb = (*clusters)[0].get_set_of_points()->size();
    for (unsigned int j = 1; j < (*clusters).size(); j++){
      if((*clusters)[j].get_set_of_points()->size() > maxim_numb){
        maxim_numb = (*clusters)[j].get_set_of_points()->size();
        biggest_index = j;
      }
    }

    //prepei na sigoureutw oti den yparxei adeio cluster. An yparxei, pare to makrynotero shmeio tou pio megalou cluster. Einai apodekth methodos
    for (unsigned int j = 0; j < (*clusters).size(); j++){
      std::unordered_map<std::string, my_vector<T> * > * clust_points = (*clusters)[j].get_set_of_points();
      if((* clust_points).size() < 1){
        //std::cout << "yparxei adeio cluster!\n";
        std::unordered_map<std::string, my_vector<T> * > * big_points = (*clusters)[biggest_index].get_set_of_points();
        double max_distn = std::numeric_limits<double>::max();           //apeiro kai kala
        max_distn = -1 * max_distn;                               //arxika -apeiro
        std::string worst_id;
        for(auto xy: *big_points){
          double dis_dist = manhattan_distance(xy.second->get_v(), (*clusters)[biggest_index].get_center_coords());
          if(dis_dist > max_distn){
            max_distn = dis_dist;
            worst_id = xy.first;
          } //telos if apostashs
        } //telos for anazhthshs sto pio xontrouli cluster
        (*clusters)[biggest_index].discorporate_point(&((*vectors_array)[worst_id]));
        objective_function -= max_distn;
        (*clusters)[j].incorporate_point(&((*vectors_array)[worst_id]));
        objective_function += manhattan_distance((*vectors_array)[worst_id].get_v(), (*clusters)[j].get_center_coords());;
      }//telos if periptwshs adeiou cluster
    } //telos for gia anazhthsh adeiwn clusters

    return objective_function;
}//telos sunarthshs




//////////////////////CURVES////////////////////////////////////
//LLOYD - each input vector to it nearest! - curve time
template <typename T>
double lloyd_ass_curve(std::vector<curve_cluster<T>>* clusters, std::unordered_map<std::string, curve<T> > *curves_array)
{
    double objective_function = 0.0;
    for (auto x :(*curves_array))
    {
        double min1 = std::numeric_limits<double>::max(); //apeiro
        int min_clust_index = -1; //to index tou cluster opou anoikei kathe x
        for (unsigned int j = 0; j < (*clusters).size(); j++)
        {
            double tmp = 0.0;
            tmp = dtw( &((*curves_array)[x.first]), (*clusters)[j].get_center_ptr()); //dist(x.second.get_v(), (*clusters)[j].get_center_coords())
            if (tmp < min1)
            {
                min1 = tmp;
                min_clust_index = j;
            }
            //std::cout << "Eimai " << x.first << " & sugkrinw " << (*clusters)[j].get_center_id() << " dist " << tmp << '\n';
        }
        (*clusters)[min_clust_index].incorporate_point(&((*curves_array)[x.first]));
        objective_function += min1;
        //std::cout << "Eimai " << x.first << " & kentro " << (*clusters)[min_clust_index].get_center_id() << " dist " << min1 << '\n';
        //std::cout << '\n';
    }
    //std::cout << "\n\n\n";
    int biggest_index = 0;
    int maxim_numb = (*clusters)[0].get_set_of_curves()->size();
    for (unsigned int j = 1; j < (*clusters).size(); j++){
      if((*clusters)[j].get_set_of_curves()->size() > maxim_numb){
        maxim_numb = (*clusters)[j].get_set_of_curves()->size();
        biggest_index = j;
      }
    }

    //prepei na sigoureutw oti den yparxei adeio cluster. An yparxei, pare to makrynotero shmeio tou pio megalou cluster. Einai apodekth methodos
    for (unsigned int j = 0; j < (*clusters).size(); j++){
      std::unordered_map<std::string, curve<T> * > * clust_points = (*clusters)[j].get_set_of_curves();
      if((* clust_points).size() < 1){
        //std::cout << "yparxei adeio cluster!\n";
        std::unordered_map<std::string, curve<T> * > * big_points = (*clusters)[biggest_index].get_set_of_curves();
        double max_distn = std::numeric_limits<double>::max();           //apeiro kai kala
        max_distn = -1 * max_distn;                               //arxika -apeiro
        std::string worst_id;
        for(auto xy: *big_points){
          double dis_dist = dtw(xy.second, (*clusters)[biggest_index].get_center_ptr());
          if(dis_dist > max_distn){
            max_distn = dis_dist;
            worst_id = xy.first;
          } //telos if apostashs
        } //telos for anazhthshs sto pio xontrouli cluster
        (*clusters)[biggest_index].discorporate_point(&((*curves_array)[worst_id]));
        objective_function -= max_distn;
        (*clusters)[j].incorporate_point(&((*curves_array)[worst_id]));
        objective_function += dtw(&((*curves_array)[worst_id]), (*clusters)[j].get_center_ptr());;
      }//telos if periptwshs adeiou cluster
    } //telos for gia anazhthsh adeiwn clusters

    return objective_function;
}



template <typename T>
my_vector<T> vectorify(curve<T> gc) //kanw to grid curve --> vector
{
    my_vector<T> gcv;
    std::vector<T> tmp;
    tmp.clear();
    for (unsigned int i = 0; i < gc.get_v_size(); i++) //kathe curve point tis gc
    {
        tmp.push_back(gc.get_points()[i].get_x());
        tmp.push_back(gc.get_points()[i].get_y());
    }
    gcv.set_v(tmp);
    if (gcv.get_v().size() != 2 * gc.get_v_size())
    {
        std::cerr << "Error pushing back grid curve's elements to vector (vectorify)\n";
        exit(-1);
    }
    gcv.set_id(gc.get_id());
    return gcv;
}


template <typename T>
void add_pad(my_vector<T> *gcv, double timi, unsigned int max_vec_size)
{
    //vale "timi" sto "gcv" mexri na exei megethos "max_vec_size"
    //std::cout << "maxdims is " << max_vec_size << std::endl;
    std::vector<T> tmp = gcv->get_v();
    //std::cout << "hello2\n";
    for (unsigned int i = 0; i < max_vec_size - gcv->get_v().size(); i++) //an px theloume max size = 10 kai exoume 2 8eloume na baloume 10-2=8 ara apo 0 ws 7
    {
        tmp.push_back(timi);
    }
    //std::cerr << "tmp size2 is " << tmp.size();
    //std::cout << "hello3\n";
    gcv->set_v(tmp);
    if (gcv->get_v().size() != max_vec_size)
    {
        std::cerr << "Error in add_pad\n";
        exit(-1);
    }
    //std::cerr << "tmp size3 is " << gcv->get_v().size();
    return; //an exei error, isws thelei na kanoume ena tmp antigrafo tou gcv, na baloume ekei ta "timi" & na epistrepsoume auto
}




template <typename T> //APO ERGASIA 1
double LSH_range_ass_curve(std::vector<curve_cluster<T>>* clusters, std::unordered_map<std::string, curve<T> > *curves_array, int number_of_grids, int number_of_lsh_hash_functions, double delta, double max_coord, double * w_done)
{
    double w = *w_done;
    if(*w_done < 0){ //den exei ypologistei, pame gia 1h kai teleutaia fora to w tou lsh
      //auto start_of_w_calc = std::chrono::high_resolution_clock::now();
      std::vector<NNpair> input_actual_NNs; //pinakas apo zeugaria actual NNs me prwto stoixeio to p
      //int n =0;
      for (auto x :(*curves_array))
      { //prepei na brw ta zeugaria ap to input gia ypologismo w

        std::string min_id1;
        double min1 = std::numeric_limits<double>::max(); //min pairnei timh apeiro
        for (auto y :(*curves_array))
        {
          if (dtw(&(x.second), &(y.second)) == 0){
            if(x.second.get_id() == y.second.get_id()){
              continue;
            }
          }

          if (dtw(&(x.second), &(y.second)) < min1)
          {
            min1 = dtw(&(x.second), &(y.second));
            min_id1 =  y.second.get_id();
          }
        }
        NNpair single_pair1(x.second.get_id(), min_id1);
        single_pair1.set_distance(min1);
        input_actual_NNs.push_back(single_pair1);
        //n++;
      }

      double tmp = 0.0;
      double mean_distance = 0;
      for (unsigned int i = 0; i < input_actual_NNs.size(); i++)
      {
        tmp += input_actual_NNs.at(i).get_distance();
      }
      //int diastaseis = vectors_array->begin().second.get_v().size(); //gia hash tables
      //std::cout << "diasteaseis " << diastaseis << "\n";
      mean_distance = tmp / input_actual_NNs.size(); //fp division
      //auto end_of_w_calc = std::chrono::high_resolution_clock::now() - start_of_w_calc;
      //long long microseconds_w = std::chrono::duration_cast<std::chrono::microseconds>(end_of_w_calc).count();
      //fprintf(stderr, "Time needed for w calculation is %lld microseconds.\n\n", microseconds_w);
      fprintf(stderr, "Value of w = %f\n", mean_distance);
      //also test gia w = 10 * mean_distance
      /*const*/ w = 4 * mean_distance; //to w pou vazw sta ai, STH XEIROTERH HARD CODED
      *w_done = w; //ypologisthke prwth kai teleutaia fora!
    }
    //std::cout << "w izz " << w << "\n";
    int n = (*curves_array).size() ;
    int Table_Size = (*curves_array).size() / 8;
    unsigned int max_dims_in = 0;

    /////////////////////////////LSH TIME////////////////////////////////
    std::unordered_map<std::string, std::pair<int, double>> owned; //to flag poy lene oi diafaneies gia to an kaparw8hke ena shmeio kai apo poio index kai me poia aktina
    std::vector<grid<T>> grids_v; //o pinakas twn grids
    for (int i = 0; i < number_of_grids; i++)
    {
      grid<T> tmp(delta, 2); //2 einai oi diastaseis mas gia tis kampyles
      grids_v.push_back(tmp);
    }//telos for gia grid init
    std::vector<my_vector<T>> input_vectors_array; //TA VECTORS POY PROEKYPSAN APO TIS KAMPYLES EISODOU
    for (int j = 0; j < number_of_grids; j++){
      input_vectors_array.clear();
      max_dims_in = 0;
      for (auto x: *curves_array){
        curve<T> grid_curve;
        grid_curve = grids_v[j].gridify(&((*curves_array)[x.first]));
        //std::cout << grid_curve.get_points().size() << std::endl;
        my_vector<T> converted_vec;
        converted_vec = vectorify(grid_curve);
        //std::cout << converted_vec.get_v().size() << std::endl;
        if (converted_vec.get_v().size() > max_dims_in)
          max_dims_in = converted_vec.get_v().size();
        input_vectors_array.push_back(converted_vec);
      }//telos for gia gridify/vectorify se ka8e kampylh
      grids_v[j].define_hash_table(Table_Size, number_of_lsh_hash_functions, max_dims_in, w); //gia auto to grid dhmiourgei ton hash table tou gia na kanoyme lsh
      for (unsigned int i = 0; i < input_vectors_array.size(); i++)
      {
        add_pad(&input_vectors_array[i], 100 * max_coord, max_dims_in);
      }


      for (unsigned int i = 0; i < input_vectors_array.size(); i++) //kanoume isomhkh ola ta vectors
      {//hasharw ta vectors k arxikopoiw flags
        //Twra kanoume LSH apo A erwthma sta vectors auta
        grids_v[j].hash_table.hash_vector(&input_vectors_array[i], &((*curves_array)[input_vectors_array[i].get_id()]));
        //analoga to ti epistrefoun apothikeuoume sta hash tables (posa?) tis antistoixes kampules
        std::pair<int, double> index_and_radius;
        index_and_radius.first = -1;
        index_and_radius.second = 0.0;
        owned[((*curves_array)[input_vectors_array[i].get_id()]).get_id()] = index_and_radius;
      }
    }//telos for gia ka8e grid


    int num_unassigned = n; //posa exoun meinei xwris anathesh
    double radius = initialize_radius_curve(clusters);
    std::vector<std::string> this_center_neighbs;
    std::vector<std::string> this_HT_neighbs;
    bool repetition = false;
    int kill_countdown = 15; //an den exoun ginei nees anatheseis meta apo tosous diplasiasmous aktinas, stop
    //std::cout << "aaaanteksaaaaa\n";
    int num_unassigned_prev = num_unassigned; //arithmos unassigned shmeiwn prin th loypa gia na sugkrinoyme proodo kathe fora kai na stamatame

    while((num_unassigned > n/10) && (kill_countdown >0) ){ //h anazhthsh range search lsh tha ginetai mexri to 90% twn shmeiwn ginei assign se kapoio kentro. Epeita klassikh methodos opws prota8hke

      num_unassigned_prev = num_unassigned;
      for(unsigned int i=0; i< clusters->size(); i++){ //gia kathe kentro twn clusters

        this_center_neighbs.clear();
        for (int j = 0; j < number_of_grids; j++){ //LSH se L hashtables
          this_HT_neighbs.clear();
          //vectorization kentrou
          curve<T> grid_center;
          grid_center = grids_v[j].gridify((*clusters)[i].get_center_ptr());
          //std::cout << grid_curve.get_points().size() << std::endl;
          my_vector<T> converted_center;
          converted_center = vectorify(grid_center);
          add_pad(&converted_center, 100 * max_coord, max_dims_in);
          this_HT_neighbs = grids_v[j].hash_table.hash_query(&converted_center,(*clusters)[i].get_center_ptr(), radius, repetition);
          //std::cout << this_HT_neighbs.size() << "-";
          this_center_neighbs.insert(this_center_neighbs.end(), this_HT_neighbs.begin(), this_HT_neighbs.end());
          //std::cout << this_center_neighbs.size() << " ";
        }
        //pros8hkh shmeiwn se cluster kai flag gia na mhn to paroyn kai ta ypoloipa clusters
        //std::cout << "eimai to cl " << (*clusters)[i].get_center_id() << "kai "<<this_center_neighbs.size() << "\n";
        for(unsigned int z=0; z< this_center_neighbs.size(); z++){
          //std::cout << owned[this_center_neighbs[z]].first ;
          //std::cout << this_center_neighbs[z] ;
          if( owned[this_center_neighbs[z]].first == -1 ){ //den exei kaparw8ei
              //std::cout << "kapakap ";
              (*clusters)[i].incorporate_point(&((*curves_array)[this_center_neighbs[z]]));
              owned[this_center_neighbs[z]].first = i; //to kaparwse
              owned[this_center_neighbs[z]].second = radius; //to kaparwse entos aktinas toshs
              num_unassigned--;
          }
          else{ //sugkrinoume me auton poy to exei kaparwsei
            if(owned[this_center_neighbs[z]].second >= radius){ //an kapoios allos to exei kaparwsei me mikroterh aktina, apofeugoume th sugkrish kai proxwrame
              if(dtw( &((*curves_array)[this_center_neighbs[z]]) , (*clusters)[i].get_center_ptr()) < dtw( &((*curves_array)[this_center_neighbs[z]]) , (*clusters)[owned[this_center_neighbs[z]].first].get_center_ptr()) ){
                (*clusters)[owned[this_center_neighbs[z]].first].discorporate_point(&((*curves_array)[this_center_neighbs[z]])); //to bgazei ap to palio
                (*clusters)[i].incorporate_point(&((*curves_array)[this_center_neighbs[z]])); //to vazei sto neo
                owned[this_center_neighbs[z]].first = i; //to kaparwse
                owned[this_center_neighbs[z]].second = radius; //to kaparwse entos aktinas toshs
              } //telos if gia sugkrish apostasewn
            } //telos if poy afora an ena kaparwmeno shmeio exei kaparw8ei apo mikroterh aktina ara den exei nohma na koitaksoume pali
          } //telos else poy afora to an ena shmeio einai kaparwmeno h oxi
        } //telos gor gia auta poy brhke auto to cluster gia authn thn aktina
      } //telos for gia ta clusters

      if(num_unassigned_prev - num_unassigned <=0 ) //den kaname proodo, arxise antistrofh metrhsh
        kill_countdown--;
      else
        kill_countdown = 15; //eixame proodo, mhdenise thn antistrofh metrhsh

      if(kill_countdown<=0) //den yphrkse veltiwsh gia sunexomenes loypes, telos
        break;

      radius = radius*2; //diplasiazoume aktina kai sunexizoume
      repetition = true;
      //std::cout << num_unassigned << "\n";
    } //telos ths while poy diplasiazoume thn aktina


    //twra an kapoio exei meinei akaparwto, prepei na paei sto kontinotero tou kentro
    double objective_function = 0.0; //gia ton elegxo metavolhs antikeimenikhs sunarthshs
    for(auto x: owned){
      if(x.second.first == -1){ //akaparwto
        double min1 = std::numeric_limits<double>::max(); //apeiro
        int min_clust_index = -1; //to index tou cluster opou anoikei kathe x
        //(*vectors_array)[x.first].get_id()
        for (unsigned int j = 0; j < (*clusters).size(); j++)
        {
            double tmp = 0.0;
            tmp = dtw( &((*curves_array)[x.first]), (*clusters)[j].get_center_ptr()); //dist(x.second.get_v(), (*clusters)[j].get_center_coords())
            if (tmp < min1)
            {
                min1 = tmp;
                min_clust_index = j;
            }
            //std::cout << "Eimai " << x.first << " & sugkrinw " << (*clusters)[j].get_center_id() << " dist " << tmp << '\n';
        }
        (*clusters)[min_clust_index].incorporate_point(&((*curves_array)[x.first]));
        objective_function += min1;
      }
      else{ //kaparwmeno, aplws ypologismos apostashs gia objective funct
        objective_function += dtw( &((*curves_array)[x.first]) , (*clusters)[x.second.first].get_center_ptr() ); //apostash kampylhs apo to kentro ths
      }
    }//telos for gia akaparwta

    int biggest_index = 0;
    int maxim_numb = (*clusters)[0].get_set_of_curves()->size();
    for (unsigned int j = 1; j < (*clusters).size(); j++){
      if((*clusters)[j].get_set_of_curves()->size() > maxim_numb){
        maxim_numb = (*clusters)[j].get_set_of_curves()->size();
        biggest_index = j;
      }
    }

    //prepei na sigoureutw oti den yparxei adeio cluster. An yparxei, pare to makrynotero shmeio tou pio megalou cluster. Einai apodekth methodos
    for (unsigned int j = 0; j < (*clusters).size(); j++){
      std::unordered_map<std::string, curve<T> * > * clust_points = (*clusters)[j].get_set_of_curves();
      if((* clust_points).size() < 1){
        //std::cout << "yparxei adeio cluster!\n";
        std::unordered_map<std::string, curve<T> * > * big_points = (*clusters)[biggest_index].get_set_of_curves();
        double max_distn = std::numeric_limits<double>::max();           //apeiro kai kala
        max_distn = -1 * max_distn;                               //arxika -apeiro
        std::string worst_id;
        for(auto xy: *big_points){
          double dis_dist = dtw(xy.second, (*clusters)[biggest_index].get_center_ptr());
          if(dis_dist > max_distn){
            max_distn = dis_dist;
            worst_id = xy.first;
          } //telos if apostashs
        } //telos for anazhthshs sto pio xontrouli cluster
        (*clusters)[biggest_index].discorporate_point(&((*curves_array)[worst_id]));
        objective_function -= max_distn;
        (*clusters)[j].incorporate_point(&((*curves_array)[worst_id]));
        objective_function += dtw(&((*curves_array)[worst_id]), (*clusters)[j].get_center_ptr());;
      }//telos if periptwshs adeiou cluster
    } //telos for gia anazhthsh adeiwn clusters

    return objective_function;
}//telos sunarthshs
